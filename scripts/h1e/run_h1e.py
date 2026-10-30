#!/usr/bin/env python3
"""H1e pairs (ANALYSIS_PLAN §6 X4): AF-AF vs PDB-PDB divergence of 1:1 orthologs with experimental structures on both
sides, plus the X4 floor (same protein, AF vs PDB). Single entry point, every sub-step guarded by a sentinel.

  python scripts/h1e/run_h1e.py --orthoxml X --root ROOT --threads N [--max-chains K] [--n-shards S]

Steps (sentinels ROOT/data/interim/h1e_run/.step<i>_<name>.done; a rerun skips finished steps):
 1 chains   Qualifying chain per UniProt accession of the OrthoXML, plan rule applied explicitly from the raw SIFTS /
            RCSB tables (data/raw/pdb_meta): rcsb_method == "X-ray", resolution <= 3.0 A, initial release date >
            2018-04-30, observed coverage = min(SIFTS-mapped residues, RCSB modeled residues) / UniProt length
            (protein_table len_2026) >= 0.70. Best chain = coverage, then resolution, then newest release
            (same ranking as scripts/h1e/build_pdb_chains.py). -> h1e_run/chains_selected.tsv.
            Why not data/processed/h1e/pdb_chains_post_af2.tsv: it ranks X-ray and EM together (so the best post-AF2
            chain may be EM while a qualifying X-ray chain exists) and its `qualifies` is null where the UniProt
            length was unknown. The agreement with that table is reported.
 2 pairs    Candidate pairs = genes with a qualifying chain, same top-level HOG, different species; labelled on
            demand with build_hogs.HogIndex (== hog_pair_labels.label_pairs); kept iff relation == ortholog and
            one_to_one == 1. --max-chains K (tests only) keeps a seeded random subset of pairs whose union of genes
            has <= K proteins. Writes h1e_run/pairs_h1e.tsv with three pair sets:
              h1e_af    gene_a/gene_b = UniProt acc (AFDB v4 models)
              h1e_pdb   gene_a/gene_b = PDB store entry "<acc>_<pdb>_<authchain>"
              h1e_afpdb (contract addition) same protein, PDB entry vs AF model; hog = NA, relation = "identical"
 3 fetch    PDB chains (scripts/h1e/fetch_pdb_chains.py, UniProt numbering, chain 'A') -> h1e_run/chains/
 4 af       AFDB v4 models: read from the local panel FoldComp DB (data/interim/structures/panel_v4 of the project),
            missing ones fetched with scripts/pairs/fetch_structures.py (keys from afdb_v4_lookup.parquet);
            decompressed to h1e_run/af/AF-<acc>-F1-model_v4.pdb
 5 store    One combined store (AF + PDB entries; strucevo.pairs.store_files) -> h1e_run/store/
 6 run      scripts/pairs/run_pairs.py --alns seqonly tdiaa, then tmalign -> ROOT/data/interim/pairs_h1e/{aln}/
 7 summary  ROOT/results/h1e/h1e_run_summary.json
Pairs whose structures are missing (no AF model, chain download/extraction failed, store rejection) stay in the
table as fail_align rows (never dropped).
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import os
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import polars as pl

PROJ = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJ / "src"))
sys.path.insert(0, str(PROJ / "scripts/h1e"))
sys.path.insert(0, str(PROJ / "scripts/hogs"))
sys.path.insert(0, str(PROJ / "scripts/ancestra"))
sys.path.insert(0, str(PROJ / "tools/ancestra-mvp-ai"))

CUT = "2018-04-30"; MAXRES = 3.0; MINCOV = 0.70
SEED = 20261030; STAGE_CAP = 16  # seeded test cap only (not in config/seeds.yaml; never used in production)
T0 = time.time()


def log(*a):
    print(f"[{time.time() - T0:7.1f}s]", *a, flush=True)


def taxid_of(leaf: str) -> int:
    return int(leaf.rsplit("_", 1)[1])


class Steps:
    def __init__(self, d: Path):
        self.d = d; self.timing = {}
        p = d / "step_timings.json"
        if p.exists():
            self.timing = json.load(open(p))

    def done(self, i, name):
        return (self.d / f".step{i}_{name}.done").exists()

    def mark(self, i, name, t, **info):
        self.timing[name] = {"s": round(t, 1), **info}
        json.dump(self.timing, open(self.d / "step_timings.json", "w"), indent=1)
        (self.d / f".step{i}_{name}.done").write_text(json.dumps(self.timing[name]))


# ------------------------------------------------------------------------------------------ step 1
def select_chains(genes: set, raw: Path, lengths: dict) -> pl.DataFrame:
    s = pl.read_csv(raw / "pdb_chain_uniprot.tsv.gz", separator="\t", skip_rows=1, infer_schema_length=0)
    s = s.rename({"PDB": "pdb_id", "CHAIN": "chain", "SP_PRIMARY": "uniprot"}).filter(
        pl.col("uniprot").is_in(list(genes))).with_columns(pl.col("SP_BEG").cast(pl.Int64), pl.col("SP_END").cast(pl.Int64))
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
    n0 = c["uniprot"].n_unique()
    rule = ((pl.col("rcsb_method") == "X-ray") & (pl.col("resolution") <= MAXRES) & (pl.col("release_date") > CUT)
            & (pl.col("coverage") >= MINCOV))
    q = (c.filter(rule).sort(["uniprot", "n_mapped", "resolution", "release_date"], descending=[False, True, False, True])
         .unique("uniprot", keep="first", maintain_order=True))
    log(f"chains: {n0} accessions (of {len(genes)}) with any SIFTS chain -> {q.height} with a qualifying chain "
        f"(X-ray, <= {MAXRES} A, released > {CUT}, coverage >= {MINCOV})")
    return q.select("uniprot", "pdb_id", "chain", "rcsb_method", "resolution", "release_date", "n_sifts",
                    "modeled_residue_count", "n_mapped", "uniprot_length", "coverage")


# ------------------------------------------------------------------------------------------ step 2
def candidate_pairs(ox, qual: set, max_chains: int | None):
    from ancestra.orthoxml import Group
    from build_hogs import HogIndex
    from hog_pair_labels import parent_map
    par = parent_map(ox.species_tree)
    sp_of = {k: v.species for k, v in ox.genes.items()}
    rows, cnt = [], Counter()
    for top in ox.groups:
        if not isinstance(top, Group):
            continue
        g = [x for x in top.iter_genes() if x in qual]
        if len(g) < 2:
            continue
        hi = HogIndex(top, sp_of, par)
        for i in range(len(g)):
            for j in range(i + 1, len(g)):
                x, y = g[i], g[j]
                if sp_of[x] == sp_of[y]:
                    continue
                cnt["candidate"] += 1
                rel, lvl, frm, o2o, dcg = hi.label(x, y)
                cnt[f"{rel}_o2o{o2o}"] += 1
                if rel == "ortholog" and o2o == 1:
                    rows.append((x, y, top.id, taxid_of(sp_of[x]), taxid_of(sp_of[y]), rel, o2o, dcg, lvl, frm))
    df = pl.DataFrame(rows, orient="row", schema=["acc_a", "acc_b", "hog", "species_a", "species_b", "relation",
                                                  "one_to_one", "dcg_id", "event_level", "dup_branch_from"])
    log(f"pairs: {dict(cnt)} -> {df.height} 1:1 ortholog pairs")
    info = {"label_counts": dict(cnt), "n_pairs_1to1": df.height}
    if max_chains and df.height:
        rng = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(STAGE_CAP,)))
        keep, genes = [], set()
        for r in rng.permutation(df.height):
            a, b = df["acc_a"][int(r)], df["acc_b"][int(r)]
            new = {a, b} - genes
            if len(genes) + len(new) <= max_chains:
                keep.append(int(r)); genes |= new
        df = df.with_row_index('_r').filter(pl.col('_r').is_in(keep)).drop('_r')
        log(f"pairs: --max-chains {max_chains} -> {df.height} pairs, {len(genes)} proteins")
        info["n_pairs_after_cap"] = df.height
    return df, info


# ------------------------------------------------------------------------------------------ step 4
def af_models(accs: list, out: Path, local_fc: Path, lookup_pq: Path, threads: int) -> dict:
    import foldcomp
    out.mkdir(parents=True, exist_ok=True)
    want = [a for a in accs if not (out / f"AF-{a}-F1-model_v4.pdb").exists()]
    lk = {}
    with open(str(local_fc) + ".lookup") as fh:
        for l in fh:
            k, nm, _ = l.split("\t"); lk[nm] = int(k)
    need = {f"AF-{a}-F1-model_v4": a for a in want}
    keys = {lk[n]: n for n in need if n in lk}
    n_local = 0
    if keys:
        recs = []
        with open(str(local_fc) + ".index") as fh:
            for l in fh:
                k, o, n = map(int, l.split("\t"))
                if k in keys:
                    recs.append((o, n, keys[k]))
        recs.sort()
        with open(local_fc, "rb") as fh:
            for o, n, nm in recs:
                fh.seek(o); _, pdb = foldcomp.decompress(fh.read(n))
                (out / f"{nm}.pdb").write_text(pdb); n_local += 1
    miss = [a for n, a in need.items() if n not in lk]
    n_remote = 0
    if miss:
        lt = (pl.scan_parquet(lookup_pq).filter(pl.col("acc").is_in(miss)).select("acc", "fc_key", "entry")
              .collect().unique("acc").sort("fc_key"))
        if lt.height:
            tsv = out.parent / "af_remote_accessions.tsv"; lt.write_csv(tsv, separator="\t")
            fc = out.parent / "af_remote_fc"
            subprocess.run([sys.executable, str(PROJ / "scripts/pairs/fetch_structures.py"), str(tsv), str(fc),
                            "--procs", "2", "--threads", str(4 * threads)], check=True)
            with open(str(fc) + ".lookup") as fl, open(str(fc) + ".index") as fi, open(fc, "rb") as fh:
                nm_of = {int(l.split("\t")[0]): l.split("\t")[1] for l in fl}
                for l in fi:
                    k, o, n = map(int, l.split("\t"))
                    fh.seek(o); _, pdb = foldcomp.decompress(fh.read(n))
                    (out / f"{nm_of[k].replace('.cif.gz', '').replace('.pdb', '')}.pdb").write_text(pdb); n_remote += 1
    have = [a for a in accs if (out / f"AF-{a}-F1-model_v4.pdb").exists()]
    return {"n_requested": len(accs), "n_local_foldcomp": n_local, "n_remote_fetched": n_remote,
            "n_available": len(have), "n_without_afdb_v4": len(accs) - len(have)}


# ------------------------------------------------------------------------------------------ step 7
def summarise(pairs: pl.DataFrame, out: Path, run_info: dict) -> dict:
    S = {"run": run_info}
    af = pairs.filter(pl.col("pair_set") == "h1e_af")
    S["pairs_by_species_pair"] = {f"{min(a, b)}-{max(a, b)}": n for (a, b), n in Counter(
        zip(af["species_a"].cast(pl.Int64).to_list(), af["species_b"].cast(pl.Int64).to_list())).items()}
    S["pairs_by_species_pair"] = dict(sorted(S["pairs_by_species_pair"].items(), key=lambda kv: -kv[1]))
    S["n_pairs_by_set"] = dict(Counter(pairs["pair_set"].to_list()))
    S["status"], S["medians"], S["paired_pdb_minus_af"] = {}, {}, {}
    for aln in ("seqonly", "tdiaa", "tmalign"):
        fs = sorted((out / aln).glob("shard_*.parquet"))
        if not fs:
            continue
        d = pl.concat([pl.read_parquet(f) for f in fs])
        S["status"][aln] = {ps: dict(Counter(g["status"].to_list())) for (ps,), g in d.group_by("pair_set")}
        ok = d.filter(pl.col("status") == "ok")
        S["medians"][aln] = {ps: {c: (None if g[c].drop_nans().len() == 0 else round(float(g[c].drop_nans().median()), 4))
                                  for c in ("aa_ml", "tdi_ml", "lddt_global", "tm_short", "n_core")} | {"n_ok": g.height}
                             for (ps,), g in ok.group_by("pair_set")}
        # paired PDB-PDB minus AF-AF (same ortholog pair, both ok)
        key = pl.col("hog") + "|" + pl.min_horizontal(pl.col("gene_a").str.split("_").list.first(),
                                                      pl.col("gene_b").str.split("_").list.first()) + "|" + \
            pl.max_horizontal(pl.col("gene_a").str.split("_").list.first(), pl.col("gene_b").str.split("_").list.first())
        a = ok.filter(pl.col("pair_set") == "h1e_af").with_columns(key.alias("k")).select("k", "aa_ml", "tdi_ml", "lddt_global")
        p = ok.filter(pl.col("pair_set") == "h1e_pdb").with_columns(key.alias("k")).select("k", "aa_ml", "tdi_ml", "lddt_global")
        j = a.join(p, on="k", suffix="_pdb")
        S["paired_pdb_minus_af"][aln] = {"n": j.height} | {
            c: (round(float((j[f"{c}_pdb"] - j[c]).drop_nans().median()), 4) if j.height else None)
            for c in ("aa_ml", "tdi_ml", "lddt_global")}
    return S


# ------------------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--orthoxml", required=True)
    ap.add_argument("--root", default=str(PROJ), help="project-layout root for outputs (default: project root)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--max-chains", type=int, default=None, help="TEST ONLY: cap on proteins (seeded pair subset)")
    ap.add_argument("--n-shards", type=int, default=None, help="default max(1, n_pairs // 300)")
    ap.add_argument("--fetch-workers", type=int, default=None, help="download threads (I/O bound; default 2*threads)")
    ap.add_argument("--no-afpdb", action="store_true", help="skip the h1e_afpdb floor set")
    ap.add_argument("--raw", default=str(PROJ / "data/raw/pdb_meta"))
    ap.add_argument("--protein-table", default=str(PROJ / "data/processed/panel/protein_table.parquet"))
    ap.add_argument("--local-foldcomp", default=str(PROJ / "data/interim/structures/panel_v4"))
    ap.add_argument("--afdb-lookup", default=str(PROJ / "data/interim/panel/afdb_v4_lookup.parquet"))
    ap.add_argument("--plan-only", action="store_true", help="steps 1-2 only (counts for the ETA)")
    a = ap.parse_args()
    root = Path(a.root).resolve()
    wd = root / "data/interim/h1e_run"; wd.mkdir(parents=True, exist_ok=True)
    pout = root / "data/interim/pairs_h1e"
    assert pout.resolve() != (root / "data/interim/pairs").resolve()
    st = Steps(wd)
    ox = None

    def load_ox():
        nonlocal ox
        if ox is None:
            from ancestra.orthoxml import read_orthoxml
            t = time.time(); ox = read_orthoxml(a.orthoxml); log(f"orthoxml: {len(ox.genes)} genes ({time.time()-t:.0f}s)")
        return ox

    # 1 chains
    if not st.done(1, "chains"):
        t = time.time(); load_ox()
        genes = set(ox.genes)
        pt = pl.read_parquet(a.protein_table, columns=["acc", "len_2026"]).filter(pl.col("acc").is_in(list(genes)))
        q = select_chains(genes, Path(a.raw), dict(zip(pt["acc"].to_list(), pt["len_2026"].to_list())))
        q.write_csv(wd / "chains_selected.tsv", separator="\t")
        info = {"n_orthoxml_genes": len(genes), "n_qualifying": q.height}
        # panel-wide (all euk180 panel proteins) qualifying count: upper bound on genes for the full-run ETA
        ptall = pl.read_parquet(a.protein_table, columns=["acc", "len_2026"])
        if ptall.height > len(genes):
            qa = select_chains(set(ptall["acc"].to_list()), Path(a.raw), dict(zip(ptall["acc"].to_list(), ptall["len_2026"].to_list())))
            info["n_qualifying_panel_all"] = qa.height
        ref = PROJ / "data/processed/h1e/pdb_chains_post_af2.tsv"
        if ref.exists():  # agreement with the H1e table (X-ray & qualifies)
            r = pl.read_csv(ref, separator="\t", infer_schema_length=0).filter(pl.col("uniprot").is_in(list(genes)))
            rq = set(r.filter((pl.col("method") == "X-ray") & (pl.col("qualifies") == "true"))["uniprot"].to_list())
            mine = set(q["uniprot"].to_list())
            info["vs_pdb_chains_post_af2"] = {"table_xray_qualifies": len(rq), "both": len(rq & mine),
                                              "only_table": len(rq - mine), "only_rule": len(mine - rq),
                                              "only_rule_table_row_is_EM": int(r.filter(pl.col("uniprot").is_in(list(mine - rq)) & (pl.col("method") == "EM")).height),
                                              "only_rule_table_qualifies_null": int(r.filter(pl.col("uniprot").is_in(list(mine - rq)) & pl.col("qualifies").is_null()).height)}
            log("agreement with pdb_chains_post_af2.tsv:", info["vs_pdb_chains_post_af2"])
        st.mark(1, "chains", time.time() - t, **info)
    q = pl.read_csv(wd / "chains_selected.tsv", separator="\t", infer_schema_length=0)
    pdbname = {u: f"{u}_{p}_{c}" for u, p, c in zip(q["uniprot"], q["pdb_id"], q["chain"])}

    # 2 pairs
    if not st.done(2, "pairs"):
        t = time.time(); load_ox()
        df, info = candidate_pairs(ox, set(pdbname), a.max_chains)
        split = {}
        sp = root / "data/interim/split.tsv"
        if sp.exists():
            s = pl.read_csv(sp, separator="\t", infer_schema_length=0)
            if {"hog", "split"} <= set(s.columns):
                split = dict(zip(s["hog"], s["split"]))
        df = df.with_columns(pl.col("hog").replace_strict(split, default="none").alias("split"))
        base = ["hog", "split", "species_a", "species_b", "relation", "one_to_one", "dcg_id", "event_level", "dup_branch_from"]
        af = df.select(pl.lit("h1e_af").alias("pair_set"), pl.col("acc_a").alias("gene_a"), pl.col("acc_b").alias("gene_b"), *base, "acc_a", "acc_b")
        pdb = df.select(pl.lit("h1e_pdb").alias("pair_set"), pl.col("acc_a").replace_strict(pdbname).alias("gene_a"),
                        pl.col("acc_b").replace_strict(pdbname).alias("gene_b"), *base, "acc_a", "acc_b")
        parts = [af, pdb]
        if not a.no_afpdb:
            g = pl.concat([df.select(pl.col("acc_a").alias("acc"), pl.col("species_a").alias("sp")),
                           df.select(pl.col("acc_b").alias("acc"), pl.col("species_b").alias("sp"))]).unique("acc").sort("acc")
            parts.append(g.select(pl.lit("h1e_afpdb").alias("pair_set"), pl.col("acc").replace_strict(pdbname).alias("gene_a"),
                                  pl.col("acc").alias("gene_b"), pl.lit("NA").alias("hog"), pl.lit("none").alias("split"),
                                  pl.col("sp").alias("species_a"), pl.col("sp").alias("species_b"),
                                  pl.lit("identical").alias("relation"), pl.lit(None, dtype=pl.Int64).alias("one_to_one"),
                                  pl.lit("").alias("dcg_id"), pl.lit("").alias("event_level"), pl.lit("").alias("dup_branch_from"),
                                  pl.col("acc").alias("acc_a"), pl.col("acc").alias("acc_b")))
        allp = pl.concat([p.with_columns(pl.col("one_to_one").cast(pl.Int64)) for p in parts])
        allp.write_csv(wd / "pairs_h1e.tsv", separator="\t")
        info["n_by_set"] = dict(Counter(allp["pair_set"].to_list()))
        info["n_genes"] = int(pl.concat([df["acc_a"], df["acc_b"]]).n_unique())
        log("pairs_h1e.tsv:", info["n_by_set"])
        st.mark(2, "pairs", time.time() - t, **info)
    pairs = pl.read_csv(wd / "pairs_h1e.tsv", separator="\t", infer_schema_length=0)
    genes = sorted(set(pairs["acc_a"].to_list()) | set(pairs["acc_b"].to_list()))
    if a.plan_only:
        log("plan-only:", json.dumps(st.timing)); return

    # 3 fetch PDB chains
    if not st.done(3, "fetch"):
        import fetch_pdb_chains as F
        t = time.time()
        cd = wd / "chains"; cd.mkdir(exist_ok=True); cache = wd / "cache"; cache.mkdir(exist_ok=True)
        qi = {u: (p, c) for u, p, c in zip(q["uniprot"], q["pdb_id"], q["chain"])}
        jobs = [(g, qi[g][0], qi[g][1], str(cd), str(cache)) for g in genes]
        res = []
        with cf.ThreadPoolExecutor(a.fetch_workers or 2 * a.threads) as ex:
            for i, r in enumerate(ex.map(F.one, jobs)):
                res.append(r)
                if (i + 1) % 100 == 0:
                    log(f"fetch {i+1}/{len(jobs)} ({(i+1)/(time.time()-t):.2f} chains/s)")
        keys = ["uniprot", "pdb_id", "chain", "n_extracted", "n_sifts_mapped", "n_match", "identity", "uniprot_length",
                "up_first", "up_last", "cov_extracted", "file", "status"]
        pl.DataFrame([{k: str(r.get(k, "")) for k in keys} for r in res]).write_csv(cd / "fetch_summary.tsv", separator="\t")
        dt = time.time() - t
        nok = sum(r["status"] == "ok" for r in res)
        fails = Counter(r["status"][:60] for r in res if r["status"] != "ok")
        log(f"fetch: {nok}/{len(res)} ok in {dt:.0f}s ({len(res)/max(dt,1e-9):.2f} chains/s); failures {dict(fails)}")
        st.mark(3, "fetch", dt, n=len(res), n_ok=nok, chains_per_s=round(len(res) / max(dt, 1e-9), 3),
                failures=dict(fails))

    # 4 AF models
    if not st.done(4, "af"):
        t = time.time()
        info = af_models(genes, wd / "af", Path(a.local_foldcomp), Path(a.afdb_lookup), a.threads)
        log("af:", info); st.mark(4, "af", time.time() - t, **info)

    # 5 store
    if not st.done(5, "store"):
        from strucevo.pairs.store_files import build_store_from_files
        t = time.time()
        fsum = pl.read_csv(wd / "chains/fetch_summary.tsv", separator="\t", infer_schema_length=0).filter(pl.col("status") == "ok")
        ents = [dict(path=str(wd / "af" / f"AF-{g}-F1-model_v4.pdb"), name=f"AF-{g}-F1-model_v4", acc=g, kind="af")
                for g in genes if (wd / "af" / f"AF-{g}-F1-model_v4.pdb").exists()]
        for u, f in zip(fsum["uniprot"], fsum["file"]):
            nm = Path(f).stem
            ents.append(dict(path=f, name=nm, acc=nm, kind="pdb"))
        info = build_store_from_files(ents, str(wd / "store"), a.threads)
        log("store:", info); st.mark(5, "store", time.time() - t, **{k: v for k, v in info.items() if k != "foldseek_version"})

    # 6 pairs pipeline
    nsh = a.n_shards or max(1, pairs.height // 300)
    for i, alns in ((6, ["seqonly", "tdiaa"]), (7, ["tmalign"])):
        name = "run_" + "_".join(alns)
        if st.done(i, name):
            continue
        t = time.time()
        cmd = [sys.executable, str(PROJ / "scripts/pairs/run_pairs.py"), "--pairs", str(wd / "pairs_h1e.tsv"),
               "--store", str(wd / "store"), "--out", str(pout), "--threads", str(a.threads),
               "--n-shards", str(nsh), "--alns", *alns]
        log("run:", " ".join(cmd))
        subprocess.run(cmd, check=True)
        dt = time.time() - t
        log(f"{name}: {pairs.height} pairs in {dt:.0f}s ({pairs.height/dt:.1f} pairs/s)")
        st.mark(i, name, dt, n_pairs=pairs.height, pairs_per_s=round(pairs.height / dt, 2))

    # 7 summary + ETA
    t = time.time()
    S = summarise(pairs, pout, {"orthoxml": a.orthoxml, "root": str(root), "threads": a.threads,
                                "max_chains": a.max_chains, "n_shards": nsh, "steps": st.timing})
    tm = st.timing
    ng = tm["pairs"].get("n_genes", len(genes)); npairs = pairs.height
    eta = {"note": "linear scaling of measured per-item rates; euk180 counts from --plan-only if available"}
    plan = PROJ / "data/interim/h1e_run/step_timings.json"
    if plan.exists() and root != PROJ:
        pt = json.load(open(plan)).get("pairs", {})
        if pt:
            eta["euk180_genes"] = pt.get("n_genes"); eta["euk180_pairs"] = sum(pt.get("n_by_set", {}).values())
    nq = tm["chains"].get("n_qualifying"); nqa = tm["chains"].get("n_qualifying_panel_all")
    if "euk180_genes" not in eta and nq and nqa:  # no euk180 counts yet: scale by qualifying genes (pairs ~ genes^2)
        f = nqa / nq
        eta.update(euk180_genes=round(ng * f), euk180_pairs=round(npairs * f * f),
                   scaling=f"panel-wide qualifying accessions {nqa} / this OrthoXML {nq}; pairs scaled by the square")
    G = eta.get("euk180_genes") or ng; P = eta.get("euk180_pairs") or npairs
    rate_fetch = tm.get("fetch", {}).get("chains_per_s")
    if rate_fetch:
        eta["fetch_s"] = round(G / rate_fetch)
    eta["af_store_s"] = round((tm.get("af", {}).get("s", 0) + tm.get("store", {}).get("s", 0)) * G / max(ng, 1))
    eta["pairs_s"] = round(sum(tm.get(k, {}).get("s", 0) for k in ("run_seqonly_tdiaa", "run_tmalign")) * P / max(npairs, 1))
    eta["total_s"] = sum(v for k, v in eta.items() if k.endswith("_s"))
    S["eta_full_run"] = eta
    rd = root / "results/h1e"; rd.mkdir(parents=True, exist_ok=True)
    json.dump(S, open(rd / "h1e_run_summary.json", "w"), indent=1, default=str)
    log("summary ->", rd / "h1e_run_summary.json"); log("ETA:", eta)
    (wd / ".step8_summary.done").write_text(json.dumps({"s": round(time.time() - t, 1)}))


if __name__ == "__main__":
    main()
