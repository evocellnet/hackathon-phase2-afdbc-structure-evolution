"""Structure store (strucevo.pairs.store.Store format) built from structure FILES: AF models and PDB chains together.

Used by the H1e experimental-structure control (ANALYSIS_PLAN §6 X4; scripts/h1e/run_h1e.py).

    build_store_from_files(entries, out_dir, threads)
    entries: list of dicts {path, name, acc, kind}  kind = "af" (pLDDT from the CA B-factor) or "pdb" (experimental)

Layout is identical to ``store.build_store`` (meta.parquet, ca/plddt/aa/tdi .npy, Foldseek ``db``), plus a ``kind``
column in meta.parquet and ``rejected.tsv``. Names must be unique; Foldseek names the entry after the file stem, so
every input is staged as ``<out>/input/<name>.pdb`` (symlink or decompressed AF model).

Residue bookkeeping (the arrays must match Foldseek's residue list exactly):
  * residues are read with gemmi from the first model / first chain, polymer amino-acid residues only, in file order;
  * one-letter codes from gemmi's residue table (non-standard -> parent, unknown -> X);
  * the entry is kept only if this string equals the Foldseek AA string and the 3Di string has the same length;
    otherwise it is rejected with a reason (``rejected.tsv``), never silently re-indexed;
  * residues without a CA atom are dropped before the comparison (Foldseek cannot place them either); if that
    does not reproduce Foldseek's string the entry is rejected.
pLDDT: AF entries carry the model's per-residue pLDDT. Experimental (PDB) residues get a pseudo-pLDDT of 100, so
the pairs pipeline's core mask (pLDDT >= 70 in both proteins) reduces to "all aligned columns" for PDB-PDB pairs and
to "AF pLDDT >= 70" for AF-PDB pairs. Missing (unobserved) PDB residues are simply absent from the entry, so a
PDB chain is a gapped sub-sequence of the UniProt sequence; Foldseek computes 3Di across chain breaks as if the
neighbours were adjacent (documented limitation; affects only the few residues flanking a gap).
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import time

import numpy as np

from .store import read_mmseqs_db


def _residues(path: str):
    import gemmi
    st = gemmi.read_structure(path)
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    ch = st[0][0]
    seq, ca, bf, noca = [], [], [], 0
    for r in ch:
        info = gemmi.find_tabulated_residue(r.name)
        aa_ = info is not None and info.is_amino_acid()
        if not aa_:
            # non-tabulated / non-standard polymer residue (e.g. 6M6): Foldseek keeps any residue with a backbone,
            # coded X -> keep it iff N, CA, C are present
            if not all(r.find_atom(x, "*") is not None for x in ("N", "CA", "C")):
                continue
        a = r.find_atom("CA", "*")
        if a is None:
            noca += 1
            continue
        c = info.one_letter_code.upper() if aa_ else "X"
        seq.append(c if c.isalpha() and c != " " else "X")
        ca.append((a.pos.x, a.pos.y, a.pos.z))
        bf.append(a.b_iso)
    return "".join(seq), np.asarray(ca, np.float32).reshape(-1, 3), np.asarray(bf, np.float32), noca


def build_store_from_files(entries: list[dict], out_dir: str, threads: int = 4) -> dict:
    """Idempotent (skips if BUILD.json exists)."""
    import polars as pl
    os.makedirs(out_dir, exist_ok=True)
    if os.path.exists(os.path.join(out_dir, "BUILD.json")):
        return json.load(open(os.path.join(out_dir, "BUILD.json")))
    t0 = time.time()
    names = [e["name"] for e in entries]
    if len(set(names)) != len(names):
        raise ValueError("store entry names must be unique")
    inp = os.path.join(out_dir, "input")
    os.makedirs(inp, exist_ok=True)
    for e in entries:
        dst = os.path.join(inp, e["name"] + ".pdb")
        if not os.path.lexists(dst):
            os.symlink(os.path.abspath(e["path"]), dst)
    db = os.path.join(out_dir, "db")
    subprocess.run(["foldseek", "createdb", inp, db, "--threads", str(threads), "-v", "1"], check=True)
    t_createdb = time.time() - t0
    lookup = {}
    with open(db + ".lookup") as fh:
        for l in fh:
            k, name, _ = l.rstrip("\n").split("\t")
            lookup[name] = int(k)
    aa = read_mmseqs_db(db)
    ss = read_mmseqs_db(db + "_ss")
    rows, rej = [], []
    ca_l, pl_l, aa_l, ss_l = [], [], [], []
    for e in entries:
        k = lookup.get(e["name"])
        if k is None:
            rej.append((e["name"], e["kind"], "not_in_foldseek_db", "")); continue
        try:
            seq, ca, bf, noca = _residues(e["path"])
        except Exception as ex:  # unreadable file
            rej.append((e["name"], e["kind"], "gemmi_error", str(ex)[:80])); continue
        s_aa, s_ss = aa[k].decode(), ss[k].decode()
        if s_aa != seq:
            why = "len_mismatch" if len(s_aa) != len(seq) else "aa_mismatch"
            nd = sum(x != y for x, y in zip(s_aa, seq))
            rej.append((e["name"], e["kind"], why, f"fs={len(s_aa)} gemmi={len(seq)} diff={nd} noCA={noca}")); continue
        if len(s_ss) != len(seq):
            rej.append((e["name"], e["kind"], "3di_len_mismatch", "")); continue
        pll = bf if e["kind"] == "af" else np.full(len(seq), 100.0, np.float32)
        rows.append(dict(key=k, name=e["name"], acc=e["acc"], kind=e["kind"], length=len(seq),
                         mean_plddt=float(pll.mean()), md5=hashlib.md5(seq.encode()).hexdigest(), n_noca=noca))
        ca_l.append(ca); pl_l.append(pll.astype(np.float32))
        aa_l.append(np.frombuffer(s_aa.encode(), np.uint8)); ss_l.append(np.frombuffer(s_ss.encode(), np.uint8))
    meta = pl.DataFrame(rows)
    lens = meta["length"].to_numpy()
    meta = meta.with_columns(pl.Series("idx", np.arange(meta.height)),
                             pl.Series("offset", np.concatenate([[0], np.cumsum(lens)[:-1]]).astype(np.int64)))
    meta = meta.select("idx", "key", "name", "acc", "length", "offset", "mean_plddt", "md5", "kind", "n_noca")
    if meta["acc"].is_duplicated().any():
        raise ValueError("store acc column must be unique (pairs are looked up by acc)")
    np.save(os.path.join(out_dir, "ca.npy"), np.concatenate(ca_l))
    np.save(os.path.join(out_dir, "plddt.npy"), np.concatenate(pl_l))
    np.save(os.path.join(out_dir, "aa.npy"), np.concatenate(aa_l))
    np.save(os.path.join(out_dir, "tdi.npy"), np.concatenate(ss_l))
    meta.write_parquet(os.path.join(out_dir, "meta.parquet"))
    pl.DataFrame(rej, schema=["name", "kind", "reason", "detail"], orient="row").write_csv(
        os.path.join(out_dir, "rejected.tsv"), separator="\t")
    from collections import Counter
    info = {"n_input": len(entries), "n_store": meta.height,
            "n_store_by_kind": dict(Counter(meta["kind"].to_list())),
            "n_rejected": len(rej), "rejected_by_reason": dict(Counter(f"{r[1]}:{r[2]}" for r in rej)),
            "n_residues": int(lens.sum()), "n_entries_with_noCA_residues": int((meta["n_noca"] > 0).sum()),
            "pdb_pseudo_plddt": 100.0, "t_createdb_s": round(t_createdb, 1), "t_total_s": round(time.time() - t0, 1),
            "foldseek_version": subprocess.run(["foldseek", "version"], capture_output=True, text=True).stdout.strip()}
    json.dump(info, open(os.path.join(out_dir, "BUILD.json"), "w"), indent=1)
    return info
