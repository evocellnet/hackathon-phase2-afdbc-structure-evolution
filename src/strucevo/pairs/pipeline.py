"""Sharded, resumable per-pair pipeline (ANALYSIS_PLAN §1.2; review C9).

Output layout (``out_dir``)::

    {aln}/shard_{k:04d}.parquet       one row per pair (schema = ANALYSIS_PLAN §1.2)
    {aln}/counts_shard_{k:04d}.npz    pair_id + (N,20,20) int32 count matrices: aa, tdi (core columns),
                                      aa_hA/aa_hB/tdi_hA/tdi_hB (e6_* pairs; zeros otherwise),
                                      aa_all (seqonly only: all aligned columns -> aa_ml_seq)
    {aln}/shard_{k:04d}.done          JSON sentinel: rows in/out, parquet + npz md5, status counts, timings
A shard is skipped on rerun iff its .done exists and the parquet md5 matches. ``seqonly`` runs first because
``aa_ml_seq``/``aa_pid_seq`` (selection variables) are copied into every alignment table.
"""
from __future__ import annotations

import hashlib
import json
import multiprocessing as mp
import os
import shutil
import time

import numpy as np
import polars as pl

from . import posthoc
from .align import FS_FIELDS, cleanup_foldseek, run_foldseek
from .ml import ml_from_counts
from .store import Store

ALIGNMENTS = ("seqonly", "tmalign", "tdiaa")
MIN_CORE = 50
MIN_COV = 0.5
META_STR = ["pair_set", "hog", "split", "relation", "dcg_id", "event_level", "dup_branch_from"]
SCHEMA_ORDER = [
    "pair_id", "pair_set", "hog", "split", "gene_a", "gene_b", "species_a", "species_b", "relation", "one_to_one",
    "dcg_id", "event_level", "dup_branch_from", "len_a", "len_b", "mean_plddt_a", "mean_plddt_b", "n_aln", "n_core",
    "cov_short", "status", "aa_p", "aa_ml", "aa_ml_se", "aa_sat", "tdi_p", "tdi_ml", "tdi_ml_se", "tdi_sat",
    "aa_ml_seq", "aa_pid_seq", "lddt_global", "lddt_core", "lddt_D", "lddt_D_se", "lddt_sat", "tm_short", "rmsd_core",
    "aa_ml_hA", "aa_ml_hA_se", "aa_ml_hB", "aa_ml_hB_se", "tdi_ml_hA", "tdi_ml_hA_se", "tdi_ml_hB", "tdi_ml_hB_se",
    "boot_sub"]


def md5_file(p):
    h = hashlib.md5()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def shard_of(keys: list[str], n_shards: int) -> np.ndarray:
    return np.array([int(hashlib.md5(k.encode()).hexdigest()[:12], 16) % n_shards for k in keys], dtype=np.int64)


def prepare_pairs(pairs: pl.DataFrame, store: Store, n_shards: int) -> pl.DataFrame:
    """Orient pairs (gene_a = shorter, ties lexicographic), attach store indices, pair_id and shard."""
    for c in META_STR:
        if c not in pairs.columns:
            pairs = pairs.with_columns(pl.lit("NA" if c in ("hog",) else ("none" if c == "split" else None),
                                              dtype=pl.Utf8).alias(c))
    for c in ("species_a", "species_b"):
        if c not in pairs.columns:
            pairs = pairs.with_columns(pl.lit(None, dtype=pl.Int64).alias(c))
    if "one_to_one" not in pairs.columns:
        pairs = pairs.with_columns(pl.lit(None, dtype=pl.Int8).alias("one_to_one"))
    pairs = pairs.with_columns(pl.col("species_a").cast(pl.Int64), pl.col("species_b").cast(pl.Int64),
                               pl.col("one_to_one").cast(pl.Utf8).replace("", None).cast(pl.Int8))
    ia = np.array([store.idx_of_acc.get(g, -1) for g in pairs["gene_a"].to_list()])
    ib = np.array([store.idx_of_acc.get(g, -1) for g in pairs["gene_b"].to_list()])
    la = np.where(ia >= 0, store.len[ia], -1)
    lb = np.where(ib >= 0, store.len[ib], -1)
    ga = np.array(pairs["gene_a"].to_list(), dtype=object)
    gb = np.array(pairs["gene_b"].to_list(), dtype=object)
    swap = (la > lb) | ((la == lb) & (ga > gb))
    sw = pl.col("_swap")
    pairs = pairs.with_columns(pl.Series("_swap", swap), pl.Series("_ia", ia), pl.Series("_ib", ib)).with_columns(
        pl.when(sw).then(pl.col("gene_b")).otherwise(pl.col("gene_a")).alias("gene_a"),
        pl.when(sw).then(pl.col("gene_a")).otherwise(pl.col("gene_b")).alias("gene_b"),
        pl.when(sw).then(pl.col("species_b")).otherwise(pl.col("species_a")).alias("species_a"),
        pl.when(sw).then(pl.col("species_a")).otherwise(pl.col("species_b")).alias("species_b"),
        pl.when(sw).then(pl.col("_ib")).otherwise(pl.col("_ia")).alias("qidx"),
        pl.when(sw).then(pl.col("_ia")).otherwise(pl.col("_ib")).alias("tidx")).drop("_swap", "_ia", "_ib")
    pairs = pairs.with_columns((pl.col("gene_a") + "__" + pl.col("gene_b")).alias("pair_id"))
    dup = pairs["pair_id"].is_duplicated()
    if dup.any():
        raise ValueError(f"{int(dup.sum())} duplicated pair_id, e.g. {pairs.filter(dup)['pair_id'][0]}")
    hk = [h if (h and h != "NA") else p for h, p in zip(pairs["hog"].to_list(), pairs["pair_id"].to_list())]
    return pairs.with_columns(pl.Series("shard", shard_of(hk, n_shards)))


def _chunks(lst, n):
    return [lst[i:i + n] for i in range(0, len(lst), n)]


class Runner:
    def __init__(self, store_dir: str, out_dir: str, threads: int = 48, alpha_aa: float = 1.0,
                 alpha_tdi: float = 1.0, device: str | None = None, tmalign_fast: int = 1, log=print):
        self.store = Store(store_dir)
        posthoc.STORE = self.store
        self.out = out_dir
        self.threads = threads
        self.alpha = {"aa": alpha_aa, "tdi": alpha_tdi}
        self.device = device
        self.tmalign_fast = tmalign_fast
        self.log = log
        self.pool = None

    def __enter__(self):
        self.pool = mp.get_context("fork").Pool(self.threads)
        return self

    def __exit__(self, *a):
        self.pool.close(); self.pool.join()

    # ----------------------------------------------------------------------------------------------
    def run_shard(self, sp: pl.DataFrame, k: int, aln: str) -> dict:
        d = os.path.join(self.out, aln)
        os.makedirs(d, exist_ok=True)
        pq = os.path.join(d, f"shard_{k:04d}.parquet")
        done = os.path.join(d, f"shard_{k:04d}.done")
        if os.path.exists(done) and os.path.exists(pq):
            info = json.load(open(done))
            if info.get("parquet_md5") == md5_file(pq) and info["rows_in"] == sp.height:
                info["skipped"] = True
                return info
        t_all = time.time()
        tm = {}
        n = sp.height
        qidx = sp["qidx"].to_numpy(); tidx = sp["tidx"].to_numpy()
        halves = np.array([s.startswith("e6_") for s in sp["pair_set"].to_list()])
        has = (qidx >= 0) & (tidx >= 0)
        res = [None] * n
        t0 = time.time()
        if aln == "seqonly":
            tasks = [(i, int(qidx[i]), int(tidx[i]), bool(halves[i])) for i in np.where(has)[0]]
            for part in self.pool.imap_unordered(posthoc.worker_seqonly, _chunks(tasks, 100)):
                for i, m in part:
                    res[i] = m
            tm["align+posthoc_s"] = time.time() - t0
        else:
            wd = os.path.join(d, f"work_{k:04d}")
            hv = np.where(has)[0]
            tsv, ftm = run_foldseek(self.store.db, self.store.key[qidx[hv]], self.store.key[tidx[hv]], wd, aln,
                                    threads=self.threads, tmalign_fast=self.tmalign_fast)
            tm.update(ftm)
            t1 = time.time()
            fs = pl.read_csv(tsv, separator="\t", has_header=False, new_columns=FS_FIELDS, quote_char=None,
                             schema_overrides={"lddtfull": pl.Utf8, "cigar": pl.Utf8, "query": pl.Utf8,
                                               "target": pl.Utf8})
            name2i = self.store.idx_of_name
            fq = np.array([name2i[x] for x in fs["query"].to_list()], dtype=np.int64)
            ft = np.array([name2i[x] for x in fs["target"].to_list()], dtype=np.int64)
            row_of = {(int(a), int(b)): i for i, (a, b) in enumerate(zip(qidx, tidx)) if a >= 0 and b >= 0}
            rows_fs = np.array([row_of.get((a, b), -1) for a, b in zip(fq, ft)])
            if (rows_fs < 0).any():
                raise RuntimeError(f"foldseek returned {(rows_fs < 0).sum()} rows for pairs not requested")
            if len(np.unique(rows_fs)) != len(rows_fs):
                raise RuntimeError("foldseek returned duplicated pairs")
            tm["parse_s"] = time.time() - t1
            self.n_fs_rows = len(rows_fs)
            qs = fs["qstart"].to_list(); ts = fs["tstart"].to_list(); cg = fs["cigar"].to_list()
            lf = fs["lddtfull"].to_list(); qtm = fs["qtmscore"].to_numpy()
            tasks = [(int(rows_fs[j]), int(fq[j]), int(ft[j]), int(qs[j]), int(ts[j]), cg[j], lf[j],
                      bool(halves[rows_fs[j]]), False) for j in range(len(rows_fs))]
            t1 = time.time()
            for part in self.pool.imap_unordered(posthoc.worker_foldseek, _chunks(tasks, 100)):
                for i, m in part:
                    res[i] = m
            tm["posthoc_s"] = time.time() - t1
            tm_short = np.full(n, np.nan)
            tm_short[rows_fs] = qtm
            cleanup_foldseek(wd)
            shutil.rmtree(wd, ignore_errors=True)
        aligned = np.array([r is not None for r in res])
        empty = posthoc._empty(True, aln == "seqonly")
        res = [r if r is not None else empty for r in res]
        # ---- ML (GPU)
        t1 = time.time()
        C = {}
        for key in ("aa", "tdi"):
            C[key] = np.stack([r[f"c_{key}"] for r in res]) if n else np.zeros((0, 20, 20), np.int64)
            for h in ("hA", "hB"):
                C[f"{key}_{h}"] = np.stack([r.get(f"c_{key}_{h}", empty[f"c_{key}_{h}"]) for r in res]) if n \
                    else np.zeros((0, 20, 20), np.int64)
        if aln == "seqonly":
            C["aa_all"] = np.stack([r["c_aa_all"] for r in res]) if n else np.zeros((0, 20, 20), np.int64)
        cols = {}
        for key in ("aa", "tdi"):
            r = ml_from_counts(C[key], key, self.alpha[key], self.device)
            cols[f"{key}_p"] = r["p"]; cols[f"{key}_ml"] = r["ml"]; cols[f"{key}_ml_se"] = r["se"]
            cols[f"{key}_sat"] = r["sat"]
            for h in ("hA", "hB"):
                rh = ml_from_counts(C[f"{key}_{h}"][halves], key, self.alpha[key], self.device)
                v = np.full(n, np.nan); s = np.full(n, np.nan)
                v[halves] = rh["ml"]; s[halves] = rh["se"]
                cols[f"{key}_ml_{h}"] = v; cols[f"{key}_ml_{h}_se"] = s
        if aln == "seqonly":
            r = ml_from_counts(C["aa_all"], "aa", self.alpha["aa"], self.device)
            cols["aa_ml_seq"] = r["ml"]; cols["aa_pid_seq"] = 1.0 - r["p"]
        tm["ml_s"] = time.time() - t1
        # ---- assemble
        st = self.store
        cols.update({
            "len_a": np.where(qidx >= 0, st.len[qidx], -1), "len_b": np.where(tidx >= 0, st.len[tidx], -1),
            "mean_plddt_a": np.where(qidx >= 0, st.meta["mean_plddt"].to_numpy()[qidx], np.nan),
            "mean_plddt_b": np.where(tidx >= 0, st.meta["mean_plddt"].to_numpy()[tidx], np.nan),
            "n_aln": np.array([r["n_aln"] for r in res], np.int64),
            "n_core": np.array([r["n_core"] for r in res], np.int64),
            "cov_short": np.array([r["cov_short"] for r in res], float),
            "lddt_global": np.array([r["lddt_global"] for r in res], float),
            "lddt_core": np.array([r["lddt_core"] for r in res], float),
            "rmsd_core": np.array([r["rmsd_core"] for r in res], float),
            "tm_short": tm_short if aln != "seqonly" else np.full(n, np.nan),
            "lddt_D": np.full(n, np.nan), "lddt_D_se": np.full(n, np.nan),
        })
        status = np.where(~aligned, "fail_align",
                          np.where(cols["n_core"] < MIN_CORE, "fail_core",
                                   np.where(cols["cov_short"] < MIN_COV, "fail_cov", "ok")))
        df = sp.select("pair_id", "pair_set", "hog", "split", "gene_a", "gene_b", "species_a", "species_b",
                       "relation", "one_to_one", "dcg_id", "event_level", "dup_branch_from")
        df = df.with_columns([pl.Series(k, v) for k, v in cols.items()] + [
            pl.Series("status", status), pl.lit(None, dtype=pl.Boolean).alias("lddt_sat"),
            pl.lit(False).alias("boot_sub")])
        if aln != "seqonly":
            sq = pl.read_parquet(os.path.join(self.out, "seqonly", f"shard_{k:04d}.parquet"),
                                 columns=["pair_id", "aa_ml_seq", "aa_pid_seq"])
            df = df.join(sq, on="pair_id", how="left")
        df = df.select(SCHEMA_ORDER)
        assert df.height == n
        tmp = pq + ".tmp"
        df.write_parquet(tmp)
        os.replace(tmp, pq)
        npz = os.path.join(d, f"counts_shard_{k:04d}.npz")
        np.savez_compressed(npz + ".tmp.npz", pair_id=np.array(sp["pair_id"].to_list()),
                            **{k2: v.astype(np.int32) for k2, v in C.items()})
        os.replace(npz + ".tmp.npz", npz)
        tm["total_s"] = time.time() - t_all
        sc = dict(zip(*np.unique(status, return_counts=True)))
        info = {"shard": k, "aln": aln, "rows_in": n, "rows_out": df.height, "n_aligned": int(aligned.sum()),
                "n_with_structures": int(has.sum()), "status": {str(a): int(b) for a, b in sc.items()},
                "parquet_md5": md5_file(pq), "counts_md5": md5_file(npz), "timings": tm,
                "alpha": self.alpha, "skipped": False}
        if aln != "seqonly":
            info["foldseek_rows"] = int(self.n_fs_rows)
            info["fail_align_with_structures"] = int((has & ~aligned).sum())
        json.dump(info, open(done, "w"), indent=1)
        return info

    def run(self, pairs: pl.DataFrame, alns=ALIGNMENTS, shards=None) -> list[dict]:
        alns = [a for a in ALIGNMENTS if a in alns]
        if alns and alns[0] != "seqonly" and not os.path.isdir(os.path.join(self.out, "seqonly")):
            raise RuntimeError("seqonly must run first (provides aa_ml_seq)")
        infos = []
        ks = sorted(pairs["shard"].unique().to_list()) if shards is None else shards
        for k in ks:
            sp = pairs.filter(pl.col("shard") == k).sort("pair_id")
            for aln in alns:
                info = self.run_shard(sp, k, aln)
                self.log(f"[shard {k:04d} {aln}] rows {info['rows_in']}->{info['rows_out']} "
                         f"status {info['status']} {'(skipped)' if info['skipped'] else info['timings']}")
                infos.append(info)
        return infos
