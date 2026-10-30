#!/usr/bin/env python3
"""Goal 2-ii item 2/4 prep: subsample within-HOG labelled pairs (both members eligible = seq_match, non-fragment,
in a non-singleton AFDB v3 cluster), at most K pairs per HOG (seed 20261030), and attach sequence-only AA
divergence (aa_ml_seq, LG+G alpha=1, aa_pid_seq) computed with the pairs pipeline's seqonly code path
(scripts/hogs/seqdist.py). Distances already in data/interim/hogs/seqdist_cache_alpha1.parquet are reused.
Output: data/interim/goal2/pairs_sample_dist.parquet
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / "hogs"))
from g2_common import OUT_I, ROOT, SEED  # noqa: E402

K = int(sys.argv[1]) if len(sys.argv) > 1 else 8
THREADS = 12
mc = pl.read_parquet(OUT_I / "hog_members_clu.parquet").select("acc", "eligible", "clu")
p = pl.read_csv(OUT_I / "hog_pairs_sample.tsv", separator="\t", infer_schema_length=0,
                columns=["gene_a", "gene_b", "species_a", "species_b", "relation", "event_level", "one_to_one", "hog"])
print("labelled pairs (thinned)", p.height, flush=True)
p = (p.join(mc.rename({"acc": "gene_a", "eligible": "el_a", "clu": "clu_a"}), on="gene_a")
      .join(mc.rename({"acc": "gene_b", "eligible": "el_b", "clu": "clu_b"}), on="gene_b")
      .filter(pl.col("el_a") & pl.col("el_b")))
rng = np.random.default_rng(SEED)
p = p.with_columns(r=pl.Series(rng.random(p.height))).sort("r")
p = p.with_columns(k=pl.col("r").rank("ordinal").over("hog")).filter(pl.col("k") <= K).drop("r", "k", "el_a", "el_b")
print("subsampled pairs", p.height, "hogs", p["hog"].n_unique(), flush=True)

cache = pl.read_parquet(ROOT / "data/interim/hogs/seqdist_cache_alpha1.parquet")
c2 = pl.concat([cache, cache.rename({"a": "b", "b": "a"}).select(cache.columns)]).unique(["a", "b"])
p = p.join(c2.rename({"a": "gene_a", "b": "gene_b", "ml": "aa_ml_seq", "pid": "aa_pid_seq", "sat": "aa_sat"}),
           on=["gene_a", "gene_b"], how="left")
todo = p.filter(pl.col("aa_ml_seq").is_null() & pl.col("aa_sat").is_null())
print("from cache", p.height - todo.height, "to compute", todo.height, flush=True)

import seqdist  # noqa: E402

accs = set(todo["gene_a"]) | set(todo["gene_b"])
seqdist.SEQ.update(seqdist.read_fasta_dir(ROOT / "data/interim/ancestra_input", accs))
pairs = list(zip(todo["gene_a"], todo["gene_b"]))
t0 = time.time()
ml, pid, sat = seq_out = seqdist.seq_distances(pairs, threads=THREADS, alpha=1.0)
print(f"computed {len(pairs)} in {time.time()-t0:.0f}s", flush=True)
new = todo.select("gene_a", "gene_b").with_columns(n_ml=pl.Series(ml), n_pid=pl.Series(pid), n_sat=pl.Series(sat))
p = (p.join(new, on=["gene_a", "gene_b"], how="left")
      .with_columns(aa_ml_seq=pl.coalesce("aa_ml_seq", "n_ml"), aa_pid_seq=pl.coalesce("aa_pid_seq", "n_pid"),
                    aa_sat=pl.coalesce("aa_sat", "n_sat")).drop("n_ml", "n_pid", "n_sat"))
p.write_parquet(OUT_I / "pairs_sample_dist.parquet")
print(p.describe(), flush=True)
