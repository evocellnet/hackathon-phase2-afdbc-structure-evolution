"""Structure store: Foldseek DB built from a FoldComp subset + ragged per-residue arrays (ANALYSIS_PLAN §1.2).

Layout of a store directory ``<store>/``::

    db, db_ss, db_ca, db_h, ... Foldseek structure DB (``foldseek createdb <foldcomp> db --input-format 5``)
    meta.parquet               one row per entry: idx, key (Foldseek DB key), name (AF-<acc>-F1-model_v4),
                               acc, length, offset (into the ragged arrays), mean_plddt, md5 (of the AA sequence)
    ca.npy        (Nres, 3) float32  CA coordinates (FoldComp-decoded)
    plddt.npy     (Nres,)   float32  per-residue pLDDT (FoldComp b-factors, 0.01 resolution)
    aa.npy        (Nres,)   uint8    AA letters (ASCII) from the Foldseek DB (identical to FoldComp residues; asserted)
    tdi.npy       (Nres,)   uint8    3Di letters (ASCII) from the Foldseek ``db_ss`` DB
    BUILD.json    provenance: command, foldseek version, counts, timings

All arrays are memory-mapped on load, so forked worker processes share them.
"""
from __future__ import annotations

import hashlib
import json
import multiprocessing as mp
import os
import subprocess
import time
from dataclasses import dataclass

import numpy as np

# heavy atoms per residue in FoldComp/PDB order (N, CA, C, O, side chain); CA is atom 1 of each residue
_NATOM = dict(A=5, R=11, N=8, D=8, C=6, Q=9, E=9, G=4, H=10, I=8, L=8, K=9, M=8, F=11, P=7, S=6, T=7, W=14, Y=12, V=7)


def read_mmseqs_db(path: str) -> dict[int, bytes]:
    """Read an MMseqs/Foldseek DB (data file + .index) into {key: bytes} (trailing '\\n\\0' stripped)."""
    out = {}
    with open(path + ".index") as fi:
        idx = [tuple(map(int, l.split("\t"))) for l in fi]
    if os.path.exists(path):
        data = np.memmap(path, dtype=np.uint8, mode="r")
        for k, o, n in idx:
            out[k] = bytes(data[o:o + n]).rstrip(b"\x00").rstrip(b"\n")
    else:  # split DB (.0, .1, ...): offsets are per file in index order
        raise NotImplementedError(f"split DB not supported: {path}")
    return out


def _decode_chunk(args):
    import foldcomp
    fc_path, recs = args
    res = []
    with open(fc_path, "rb") as fh:
        for name, off, ln in recs:
            fh.seek(off)
            d = foldcomp.get_data(fh.read(ln))
            seq = d["residues"]
            nat = np.fromiter((_NATOM.get(c, 0) for c in seq), dtype=np.int64, count=len(seq))
            X = np.asarray(d["coordinates"], dtype=np.float32)
            st = np.concatenate([[0], np.cumsum(nat)[:-1]])
            if nat.sum() not in (len(X), len(X) - 1) or (nat == 0).any():
                res.append((name, seq, None, None))
                continue
            res.append((name, seq, X[st + 1], np.asarray(d["b_factors"], dtype=np.float32)))
    return res


def build_store(foldcomp_db: str, out_dir: str, threads: int = 48) -> dict:
    """FoldComp DB -> Foldseek DB + ragged arrays. Idempotent (skips if BUILD.json exists)."""
    os.makedirs(out_dir, exist_ok=True)
    if os.path.exists(os.path.join(out_dir, "BUILD.json")):
        return json.load(open(os.path.join(out_dir, "BUILD.json")))
    t0 = time.time()
    db = os.path.join(out_dir, "db")
    subprocess.run(["foldseek", "createdb", foldcomp_db, db, "--input-format", "5", "--threads", str(threads),
                    "-v", "1"], check=True)
    t_createdb = time.time() - t0
    lookup = {}
    with open(db + ".lookup") as fh:
        for l in fh:
            k, name, _ = l.rstrip("\n").split("\t")
            lookup[name] = int(k)
    aa = read_mmseqs_db(db)
    ss = read_mmseqs_db(db + "_ss")
    # FoldComp offsets
    fc_name = {}
    with open(foldcomp_db + ".lookup") as fh:
        for l in fh:
            k, name, _ = l.rstrip("\n").split("\t")
            fc_name[int(k)] = name
    recs = []
    with open(foldcomp_db + ".index") as fh:
        for l in fh:
            k, o, n = map(int, l.split("\t"))
            recs.append((fc_name[k], o, n))
    chunks = [(foldcomp_db, recs[i:i + 2000]) for i in range(0, len(recs), 2000)]
    t1 = time.time()
    with mp.get_context("fork").Pool(threads) as pool:
        dec = [r for part in pool.imap(_decode_chunk, chunks) for r in part]
    t_decode = time.time() - t1
    names, keys, lens, md5s, mpl, bad = [], [], [], [], [], []
    ca_l, pl_l, aa_l, ss_l = [], [], [], []
    for name, seq, ca, pll in dec:
        k = lookup.get(name)
        if ca is None or k is None:
            bad.append(name)
            continue
        s_aa = aa[k].decode()
        s_ss = ss[k].decode()
        if not (len(s_aa) == len(s_ss) == len(ca) == len(pll)) or s_aa != seq:
            bad.append(name)
            continue
        names.append(name); keys.append(k); lens.append(len(seq))
        md5s.append(hashlib.md5(seq.encode()).hexdigest()); mpl.append(float(pll.mean()))
        ca_l.append(ca); pl_l.append(pll)
        aa_l.append(np.frombuffer(s_aa.encode(), np.uint8)); ss_l.append(np.frombuffer(s_ss.encode(), np.uint8))
    lens = np.array(lens, dtype=np.int64)
    off = np.concatenate([[0], np.cumsum(lens)[:-1]])
    np.save(os.path.join(out_dir, "ca.npy"), np.concatenate(ca_l))
    np.save(os.path.join(out_dir, "plddt.npy"), np.concatenate(pl_l))
    np.save(os.path.join(out_dir, "aa.npy"), np.concatenate(aa_l))
    np.save(os.path.join(out_dir, "tdi.npy"), np.concatenate(ss_l))
    import polars as pl
    acc = [n.split("-")[1] if n.startswith("AF-") else n for n in names]
    pl.DataFrame({"idx": np.arange(len(names)), "key": keys, "name": names, "acc": acc, "length": lens,
                  "offset": off, "mean_plddt": mpl, "md5": md5s}).write_parquet(os.path.join(out_dir, "meta.parquet"))
    info = {"foldcomp_db": os.path.abspath(foldcomp_db), "n_foldcomp": len(recs), "n_store": len(names),
            "n_bad": len(bad), "bad_examples": bad[:20], "n_residues": int(lens.sum()),
            "t_createdb_s": round(t_createdb, 1), "t_decode_s": round(t_decode, 1),
            "t_total_s": round(time.time() - t0, 1),
            "foldseek_version": subprocess.run(["foldseek", "version"], capture_output=True, text=True).stdout.strip()}
    json.dump(info, open(os.path.join(out_dir, "BUILD.json"), "w"), indent=1)
    return info


@dataclass
class Store:
    path: str

    def __post_init__(self):
        import polars as pl
        self.meta = pl.read_parquet(os.path.join(self.path, "meta.parquet"))
        self.ca = np.load(os.path.join(self.path, "ca.npy"), mmap_mode="r")
        self.plddt = np.load(os.path.join(self.path, "plddt.npy"), mmap_mode="r")
        self.aa = np.load(os.path.join(self.path, "aa.npy"), mmap_mode="r")
        self.tdi = np.load(os.path.join(self.path, "tdi.npy"), mmap_mode="r")
        self.off = self.meta["offset"].to_numpy()
        self.len = self.meta["length"].to_numpy()
        self.key = self.meta["key"].to_numpy()
        self.name = self.meta["name"].to_list()
        self.idx_of_acc = dict(zip(self.meta["acc"].to_list(), range(self.meta.height)))
        self.idx_of_name = dict(zip(self.name, range(self.meta.height)))
        self.db = os.path.join(self.path, "db")

    def sl(self, i):
        return slice(int(self.off[i]), int(self.off[i] + self.len[i]))

    def get(self, i):
        s = self.sl(i)
        return self.ca[s], self.plddt[s], self.aa[s], self.tdi[s]

    def seq(self, i) -> str:
        return bytes(self.aa[self.sl(i)]).decode()

    def tdi_str(self, i) -> str:
        return bytes(self.tdi[self.sl(i)]).decode()
