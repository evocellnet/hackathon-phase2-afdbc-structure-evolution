#!/usr/bin/env python3
"""Extract (acc, md5, nres) for all AFDB v4 models from the TED per-chain segmentation table.

TED (Lau et al. 2024, Zenodo 13908086) lists every AFDB v4 model with an md5 of its sequence; we use it to
check, without fetching structures, whether the CURRENT UniProt canonical sequence equals the sequence that
was modelled in AFDB v4 (UniProt 2021_04).  Output: data/interim/panel/afdb_v4_ted_md5.parquet
"""
import os, time, polars as pl
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
t = time.time()
src = os.path.join(ROOT, "data/raw/ted/ted_214m_per_chain_segmentation.tsv.gz")
df = (pl.scan_csv(src, separator="\t", has_header=False, quote_char=None,
                  new_columns=["entry", "md5", "nres"] + [f"c{i}" for i in range(3, 10)],
                  schema_overrides={"entry": pl.Utf8, "md5": pl.Utf8, "nres": pl.Int32})
      .select(pl.col("entry").str.extract(r"^AF-(.+)-F1-model_v4$", 1).alias("acc"),
              pl.col("md5"), pl.col("nres"))
      .collect())
df.write_parquet(os.path.join(ROOT, "data/interim/panel/afdb_v4_ted_md5.parquet"))
print(df.height, df["acc"].null_count(), f"{time.time()-t:.0f}s")
