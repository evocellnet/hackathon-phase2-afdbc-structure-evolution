#!/usr/bin/env python3
"""Pre-compute two helper tables from raw AFDB v4 files (used by species selection and mapping).

1. data/interim/panel/afdb_v4_lookup.parquet : FoldComp afdb_uniprot_v4.lookup parsed into
   (fc_key:int64, entry:str, acc:str, frag:int)  where entry = AF-<acc>-F<frag>-model_v4.
2. data/interim/panel/afdb_v4_taxid_counts.parquet : number of AFDB v4 entries per NCBI taxId
   (from AFDB-clusters v3 file 5 = all 214M AFDB v4 entries).
"""
import os, sys, time
import polars as pl
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
out = os.path.join(ROOT, "data/interim/panel"); os.makedirs(out, exist_ok=True)
t = time.time()
lk = pl.read_csv(os.path.join(ROOT, "data/raw/foldcomp/afdb_uniprot_v4.lookup"), separator="\t",
                 has_header=False, new_columns=["fc_key", "entry", "x"], schema_overrides={"fc_key": pl.Int64, "entry": pl.Utf8, "x": pl.Int64})
lk = lk.with_columns(pl.col("entry").str.extract(r"^AF-(.+)-F(\d+)-model_v4$", 1).alias("acc"),
                     pl.col("entry").str.extract(r"^AF-(.+)-F(\d+)-model_v4$", 2).cast(pl.Int32).alias("frag"))
print("lookup rows", lk.height, "unparsed", lk["acc"].null_count(), "frag>1", (lk["frag"] > 1).sum(), f"{time.time()-t:.0f}s", flush=True)
lk.select("fc_key", "entry", "acc", "frag", "x").write_parquet(os.path.join(out, "afdb_v4_lookup.parquet"))
del lk
t = time.time()
f5 = os.path.join(ROOT, "data/raw/afdb_clusters/v3/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz")
c = (pl.scan_csv(f5, separator="\t", has_header=False, new_columns=["rep", "entry", "flag", "taxid"])
       .group_by("taxid").agg(pl.len().alias("n_afdb_v4")).collect())
c.write_parquet(os.path.join(out, "afdb_v4_taxid_counts.parquet"))
print("taxid counts", c.height, f"{time.time()-t:.0f}s")
