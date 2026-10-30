#!/usr/bin/env python3
"""X4 (ANALYSIS_PLAN §6, DESIGN §11.10) input builder: residue-matched AF-AF vs PDB-PDB ortholog pairs.

  python scripts/h1e/prep_x4.py [--threads 12] [--fetch-workers 8]

Why a second run (the node-1 run, scripts/h1e/run_h1e.py, compared unmatched structures):
  PDB chains contain only observed residues; AF models are full length. The node-1 pairs therefore compared
  different residue sets (PDB-PDB: observed residues, pseudo-pLDDT 100 -> all aligned columns; AF-AF: all residues,
  pLDDT >= 70 core), and unobserved PDB residues create chain breaks that perturb 3Di.

Matching rule (per protein x chosen PDB chain):
  M = UniProt positions whose PDB residue has N, CA and C atoms, exists in the AFDB v4 model and has the SAME amino
  acid as the AF model (engineered mutations / sequence conflicts dropped).
  * "X_<acc>_<pdb>_<ch>" = PDB coordinates of M, B-factor := AF pLDDT of the same position
  * "A_<acc>_<pdb>_<ch>" = AF coordinates of M (same residues, same order, same pLDDT)
  Both enter the store as kind "af" (pLDDT read from the B-factor), so the pipeline core mask (pLDDT >= 70 in both
  proteins) selects the same residues in AF and PDB, sequences are identical, chain breaks are identical, and
  TM-score / LDDT normalisations use the same lengths. Unmatched originals are also stored for reference:
  "F_<acc>" (full AF model) and "R_<acc>_<pdb>_<ch>" (raw PDB chain, pseudo-pLDDT 100).

Chain choices (X-ray, resolution <= 3.0 A, observed coverage >= 0.70, ranked coverage > resolution > newest):
  post = released after 2018-04-30 (DESIGN §11.10; identical rule to run_h1e.py step 1)
  any  = no date restriction
  pre  = released on or before 2018-04-30 (inside the AF2 training window)
Pairs: Ancestra euk180 1:1 orthologs (HogIndex), different species, both proteins with a chain of that choice.

Pair sets (pair_set column; one row per pair per alignment in OUT/pairs/{aln}):
  x4m_af_<c>, x4m_pdb_<c>   matched AF-AF / PDB-PDB for chain choice c
  x4m_floor_<c>             same protein: A_ vs X_ (matched AF vs own PDB)
  x4u_af, x4u_pdb_<c>       unmatched reference: F_ vs F_, R_ vs R_
Output: data/interim/h1e_x4/{chains,af,matched,store,pairs}, pairs_x4.tsv, pair_meta.tsv, match_summary.tsv.
ML distances at the tuned alpha (config/tuning.yaml: aa 1.64, tdi 1.60).
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import polars as pl
import yaml

PROJ = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJ / "scripts/h1e"))
import run_h1e as R  # noqa: E402  (select_chains, candidate_pairs, af_models; sets sys.path)

ORTHOXML = PROJ / "data/interim/ancestra_euk180/ancestra/hogs.orthoxml"
WD = PROJ / "data/interim/h1e_x4"
CUT = "2018-04-30"
AA3 = None


def log(*a):
    R.log(*a)


def choose_chains(genes, lengths):
    """Return {choice: DataFrame(uniprot, pdb_id, chain, resolution, release_date, coverage)}."""
    out = {}
    for c, cut in (("post", CUT), ("any", "0000-00-00")):
        R.CUT = cut
        out[c] = R.select_chains(genes, PROJ / "data/raw/pdb_meta", lengths)
    m = pl.read_csv(PROJ / "data/raw/pdb_meta/pdb_entry_meta.tsv.gz", separator="\t", infer_schema_length=0)
    pre_ids = set(m.filter(pl.col("initial_release_date") <= CUT)["pdb_id"].to_list())
    out["pre"] = _select_pre(genes, lengths, pre_ids)
    R.CUT = CUT
    return out


def _select_pre(genes, lengths, pre_ids):
    # copy run_h1e.select_chains logic with an extra entry filter (keeps the rule identical)
    raw = PROJ / "data/raw/pdb_meta"
    s = pl.read_csv(raw / "pdb_chain_uniprot.tsv.gz", separator="\t", skip_rows=1, infer_schema_length=0)
    s = s.rename({"PDB": "pdb_id", "CHAIN": "chain", "SP_PRIMARY": "uniprot"}).filter(
        pl.col("uniprot").is_in(list(genes)) & pl.col("pdb_id").is_in(list(pre_ids))).with_columns(
        pl.col("SP_BEG").cast(pl.Int64), pl.col("SP_END").cast(pl.Int64))
    k = ["pdb_id", "chain", "uniprot"]
    s = s.sort(k + ["SP_BEG"]).with_columns(pl.col("SP_END").cum_max().over(k).shift(1).over(k).fill_null(0).alias("prev_end"))
    s = s.with_columns(pl.max_horizontal(pl.lit(0), pl.col("SP_END") - pl.max_horizontal(pl.col("SP_BEG") - 1, pl.col("prev_end"))).alias("new_res"))
    c = s.group_by(k).agg(pl.col("new_res").sum().alias("n_sifts"))
    m = pl.read_csv(raw / "pdb_entry_meta.tsv.gz", separator="\t", infer_schema_length=0).with_columns(
        pl.col("resolution").cast(pl.Float64, strict=False))
    c = c.join(m.select("pdb_id", "rcsb_method", "resolution", pl.col("initial_release_date").alias("release_date")), on="pdb_id")
    md = pl.read_csv(raw / "pdb_chain_modeled.tsv.gz", separator="\t", schema_overrides={
        "modeled_residue_count": pl.Int64, "entity_seq_length": pl.Int64, "chain": pl.String}).unique(["pdb_id", "chain"])
    c = c.join(md.select("pdb_id", "chain", "modeled_residue_count"), on=["pdb_id", "chain"], how="left")
    c = c.with_columns(pl.min_horizontal("n_sifts", "modeled_residue_count").fill_null(pl.col("n_sifts")).alias("n_mapped"))
    c = c.join(pl.DataFrame({"uniprot": list(lengths), "uniprot_length": list(lengths.values())}), on="uniprot", how="left")
    c = c.with_columns((pl.col("n_mapped") / pl.col("uniprot_length")).alias("coverage"))
    rule = ((pl.col("rcsb_method") == "X-ray") & (pl.col("resolution") <= R.MAXRES) & (pl.col("release_date") <= CUT)
            & (pl.col("coverage") >= R.MINCOV))
    q = (c.filter(rule).sort(["uniprot", "n_mapped", "resolution", "release_date"], descending=[False, True, False, True])
         .unique("uniprot", keep="first", maintain_order=True))
    log(f"chains pre: {q.height} accessions with a qualifying pre-cutoff X-ray chain")
    return q.select("uniprot", "pdb_id", "chain", "rcsb_method", "resolution", "release_date", "n_sifts",
                    "modeled_residue_count", "n_mapped", "uniprot_length", "coverage")


def fetch_group(job):
    """All chains of one PDB entry in one thread (avoids the shared-cache race of the node-1 run)."""
    import fetch_pdb_chains as F
    return [F.one(j) for j in job]


# ------------------------------------------------------------------------------------------ residue matching
def read_res(path):
    import gemmi
    st = gemmi.read_structure(str(path)); st.remove_hydrogens(); st.remove_alternative_conformations()
    out = {}
    for r in st[0][0]:
        info = gemmi.find_tabulated_residue(r.name)
        if info is None or not info.is_amino_acid():
            continue
        if not all(r.find_atom(x, "*") is not None for x in ("N", "CA", "C")):
            continue
        out[r.seqid.num] = (info.one_letter_code.upper(), r)
    return out


def write_res(path, residues, bfac):
    import gemmi
    st = gemmi.Structure(); st.name = path.stem
    mo = gemmi.Model("1"); ch = gemmi.Chain("A")
    for (num, r), b in zip(residues, bfac):
        nr = gemmi.Residue(); nr.name = r.name; nr.seqid = gemmi.SeqId(num, " "); nr.het_flag = "A"
        for a in r:
            a2 = a.clone(); a2.b_iso = float(b); a2.occ = 1.0; a2.altloc = "\0"; nr.add_atom(a2)
        ch.add_residue(nr)
    mo.add_chain(ch); st.add_model(mo); st.setup_entities()
    st.write_pdb(str(path))


def match_one(job):
    acc, pdbf, aff, outdir = job
    tag = Path(pdbf).stem  # <acc>_<pdb>_<chain>
    try:
        P = read_res(pdbf); A = read_res(aff)
        nums = sorted(P)
        keep = [u for u in nums if u in A and A[u][0] == P[u][0]]
        n_mis = sum(1 for u in nums if u in A and A[u][0] != P[u][0])
        n_out = sum(1 for u in nums if u not in A)
        plddt = [A[u][1].find_atom("CA", "*").b_iso for u in keep]
        if len(keep) >= 30:
            write_res(Path(outdir) / f"X_{tag}.pdb", [(u, P[u][1]) for u in keep], plddt)
            write_res(Path(outdir) / f"A_{tag}.pdb", [(u, A[u][1]) for u in keep], plddt)
        return dict(acc=acc, tag=tag, n_pdb_bb=len(nums), n_af=len(A), n_keep=len(keep), n_mismatch=n_mis,
                    n_outside_af=n_out, frac_keep_of_af=round(len(keep) / len(A), 4),
                    mean_plddt_keep=round(float(np.mean(plddt)), 2) if plddt else None,
                    frac_core_keep=round(float(np.mean(np.array(plddt) >= 70)), 4) if plddt else None,
                    status="ok" if len(keep) >= 30 else "too_few")
    except Exception as e:  # noqa: BLE001
        return dict(acc=acc, tag=tag, status=f"FAIL: {e}"[:200])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--threads", type=int, default=12)
    ap.add_argument("--fetch-workers", type=int, default=8)
    a = ap.parse_args()
    WD.mkdir(parents=True, exist_ok=True)
    tun = yaml.safe_load(open(PROJ / "config/tuning.yaml"))
    alpha_aa, alpha_tdi = float(tun["alpha"]["aa"]), float(tun["alpha"]["tdi"])

    # 1 chains + pairs ------------------------------------------------------------------------------------
    if not (WD / "pairs_base.tsv").exists():
        from ancestra.orthoxml import read_orthoxml
        t = time.time(); ox = read_orthoxml(str(ORTHOXML)); log(f"orthoxml {len(ox.genes)} genes {time.time()-t:.0f}s")
        genes = set(ox.genes)
        pt = pl.read_parquet(PROJ / "data/processed/panel/protein_table.parquet", columns=["acc", "len_2026"]).filter(
            pl.col("acc").is_in(list(genes)))
        lengths = dict(zip(pt["acc"].to_list(), pt["len_2026"].to_list()))
        ch = choose_chains(genes, lengths)
        rows = []
        for c, q in ch.items():
            rows.append(q.with_columns(pl.lit(c).alias("choice")))
        chains = pl.concat(rows)
        chains.write_csv(WD / "chains_selected.tsv", separator="\t")
        union = set(chains["uniprot"].to_list())
        df, info = R.candidate_pairs(ox, union, None)
        sp = pl.read_csv(PROJ / "data/interim/split.tsv", separator="\t", infer_schema_length=0)
        split = dict(zip(sp["hog"], sp["split"]))
        df = df.with_columns(pl.col("hog").replace_strict(split, default="none").alias("split"))
        df.write_csv(WD / "pairs_base.tsv", separator="\t")
        json.dump(info, open(WD / "pairs_base_info.json", "w"), indent=1)
    chains = pl.read_csv(WD / "chains_selected.tsv", separator="\t", infer_schema_length=0)
    base = pl.read_csv(WD / "pairs_base.tsv", separator="\t", infer_schema_length=0)
    by_choice = {c: dict(zip(g["uniprot"], zip(g["pdb_id"], g["chain"]))) for (c,), g in chains.group_by("choice")}
    log("pairs_base:", base.height, {c: len(v) for c, v in by_choice.items()})

    # 2 fetch chains needed ----------------------------------------------------------------------------------
    need = set()
    for c, m in by_choice.items():
        for x, y in zip(base["acc_a"], base["acc_b"]):
            if x in m and y in m:
                need.add((x,) + m[x]); need.add((y,) + m[y])
    cdir = WD / "chains"; cdir.mkdir(exist_ok=True)
    cache = PROJ / "data/interim/h1e_run/cache"
    old = PROJ / "data/interim/h1e_run/chains"
    todo = defaultdict(list)
    for acc, pdb, chn in sorted(need):
        f = cdir / f"{acc}_{pdb}_{chn}.pdb"
        if f.exists():
            continue
        if (old / f.name).exists():
            f.symlink_to(old / f.name); continue
        todo[pdb].append((acc, pdb, chn, str(cdir), str(cache)))
    if todo:
        t = time.time(); res = []
        with cf.ThreadPoolExecutor(a.fetch_workers) as ex:
            for r in ex.map(fetch_group, list(todo.values())):
                res.extend(r)
        log(f"fetched {sum(r['status']=='ok' for r in res)}/{len(res)} chains in {time.time()-t:.0f}s;",
            dict(Counter(r["status"][:50] for r in res if r["status"] != "ok")))
        prev = pl.read_csv(cdir / "fetch_summary.tsv", separator="\t", infer_schema_length=0) if (cdir / "fetch_summary.tsv").exists() else None
        new = pl.DataFrame([{k: str(r.get(k, "")) for k in ("uniprot", "pdb_id", "chain", "identity", "status")} for r in res])
        (pl.concat([prev, new]) if prev is not None else new).write_csv(cdir / "fetch_summary.tsv", separator="\t")

    # 3 AF models --------------------------------------------------------------------------------------------
    accs = sorted({x for x, _, _ in need})
    afd = WD / "af"; afd.mkdir(exist_ok=True)
    for x in accs:
        o = PROJ / f"data/interim/h1e_run/af/AF-{x}-F1-model_v4.pdb"
        if o.exists() and not (afd / o.name).exists():
            (afd / o.name).symlink_to(o)
    info = R.af_models(accs, afd, PROJ / "data/interim/structures/panel_v4",
                       PROJ / "data/interim/panel/afdb_v4_lookup.parquet", a.threads)
    log("af:", info)

    # 4 residue matching -------------------------------------------------------------------------------------
    md = WD / "matched"; md.mkdir(exist_ok=True)
    if not (WD / "match_summary.tsv").exists():
        jobs = [(acc, str(cdir / f"{acc}_{p}_{c}.pdb"), str(afd / f"AF-{acc}-F1-model_v4.pdb"), str(md))
                for acc, p, c in sorted(need) if (cdir / f"{acc}_{p}_{c}.pdb").exists() and (afd / f"AF-{acc}-F1-model_v4.pdb").exists()]
        with cf.ProcessPoolExecutor(a.threads) as ex:
            ms = list(ex.map(match_one, jobs, chunksize=8))
        pl.DataFrame(ms).write_csv(WD / "match_summary.tsv", separator="\t")
    ms = pl.read_csv(WD / "match_summary.tsv", separator="\t", infer_schema_length=0)
    okm = set(ms.filter(pl.col("status") == "ok")["tag"].to_list())
    log("matched:", dict(Counter(ms["status"].to_list())))

    # 5 store --------------------------------------------------------------------------------------------------
    from strucevo.pairs.store_files import build_store_from_files
    ents, seen = [], set()
    for tag in sorted(okm):
        acc = tag.split("_")[0]
        ents += [dict(path=str(md / f"X_{tag}.pdb"), name=f"X_{tag}", acc=f"X_{tag}", kind="af"),
                 dict(path=str(md / f"A_{tag}.pdb"), name=f"A_{tag}", acc=f"A_{tag}", kind="af"),
                 dict(path=str(cdir / f"{tag}.pdb"), name=f"R_{tag}", acc=f"R_{tag}", kind="pdb")]
        if acc not in seen:
            seen.add(acc); ents.append(dict(path=str(afd / f"AF-{acc}-F1-model_v4.pdb"), name=f"F_{acc}", acc=f"F_{acc}", kind="af"))
    sinfo = build_store_from_files(ents, str(WD / "store"), a.threads)
    log("store:", {k: v for k, v in sinfo.items() if k != "foldseek_version"})
    from strucevo.pairs.store import Store
    instore = set(Store(str(WD / "store")).idx_of_acc)

    # 6 pair table --------------------------------------------------------------------------------------------
    rows, meta = [], []
    keep_meta = ["hog", "split", "species_a", "species_b", "relation", "one_to_one", "dcg_id", "event_level", "dup_branch_from"]
    br = base.to_dicts()
    seen_pairs = set()

    def add(ps, ga, gb, r, choice, acc_a, acc_b, kind):
        if ga not in instore or gb not in instore:
            return
        key = (ps, *sorted((ga, gb)))
        if key in seen_pairs:
            return
        seen_pairs.add(key)
        d = {"pair_set": ps, "gene_a": ga, "gene_b": gb} | {k: (r[k] if r else ("NA" if k == "hog" else None)) for k in keep_meta}
        if not r:
            d["split"] = "none"; d["relation"] = "identical"
        rows.append(d)
        meta.append(dict(pair_set=ps, g1=min(ga, gb), g2=max(ga, gb), choice=choice, acc_a=acc_a, acc_b=acc_b, kind=kind,
                         hog=(r["hog"] if r else "NA")))

    for c, m in by_choice.items():
        for r in br:
            x, y = r["acc_a"], r["acc_b"]
            if x in m and y in m:
                tx, ty = f"{x}_{m[x][0]}_{m[x][1]}", f"{y}_{m[y][0]}_{m[y][1]}"
                if tx in okm and ty in okm:
                    add(f"x4m_af_{c}", f"A_{tx}", f"A_{ty}", r, c, x, y, "af")
                    add(f"x4m_pdb_{c}", f"X_{tx}", f"X_{ty}", r, c, x, y, "pdb")
                    add(f"x4u_pdb_{c}", f"R_{tx}", f"R_{ty}", r, c, x, y, "pdb_raw")
                    add("x4u_af", f"F_{x}", f"F_{y}", r, "na", x, y, "af_full")
                    for acc, tg in ((x, tx), (y, ty)):
                        add(f"x4m_floor_{c}", f"A_{tg}", f"X_{tg}", None, c, acc, acc, "floor")
    pt = pl.DataFrame(rows).with_columns(pl.all().cast(pl.Utf8))
    log("pair sets (before de-duplication):", dict(Counter(pt["pair_set"].to_list())))
    # the same structure pair can belong to several chain choices (same chain chosen): run it once; set membership
    # is kept in pair_meta.tsv (pair_set, g1, g2) and re-expanded by analyse_h1e.py
    pt = pt.with_columns(pl.min_horizontal("gene_a", "gene_b").alias("_g1"), pl.max_horizontal("gene_a", "gene_b").alias("_g2"))
    pt = pt.unique(["_g1", "_g2"], keep="first", maintain_order=True).drop("_g1", "_g2")
    pt.write_csv(WD / "pairs_x4.tsv", separator="\t")
    pl.DataFrame(meta).write_csv(WD / "pair_meta.tsv", separator="\t")
    log("unique structure pairs:", pt.height)

    # 7 pipeline ----------------------------------------------------------------------------------------------
    for alns in (["seqonly", "tdiaa"], ["tmalign"]):
        cmd = [sys.executable, str(PROJ / "scripts/pairs/run_pairs.py"), "--pairs", str(WD / "pairs_x4.tsv"),
               "--store", str(WD / "store"), "--out", str(WD / "pairs"), "--threads", str(a.threads),
               "--n-shards", "4", "--alns", *alns, "--alpha-aa", str(alpha_aa), "--alpha-tdi", str(alpha_tdi)]
        log("run:", " ".join(cmd)); subprocess.run(cmd, check=True)
    log("done")


if __name__ == "__main__":
    main()
