#!/usr/bin/env python3
"""Overlap of AFDB entries / cluster reps between AFDB Clusters v3 and v6 releases."""
import json
from pathlib import Path
import polars as pl
ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/afdb_clusters"
cols = ["repId", "entryId", "cluFlag", "taxId"]
v3 = pl.read_csv(RAW / "v3/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz", separator="\t", has_header=False, new_columns=cols, infer_schema_length=0)
v6 = pl.read_csv(RAW / "v6/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz", separator="\t", has_header=False, new_columns=cols, infer_schema_length=0)
c3 = pl.read_csv(RAW / "v3/2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv.gz", separator="\t", has_header=False, infer_schema_length=0).select(pl.col("column_0").alias("repId"))
c6 = pl.read_csv(RAW / "v6/2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv.gz", separator="\t", has_header=False, infer_schema_length=0).select(pl.col("column_0").alias("repId"))
res = {}
e3 = v3.select("entryId").unique(); e6 = v6.select("entryId").unique()
res["v3_entries"] = e3.height; res["v6_entries"] = e6.height
res["v3_entries_in_v6"] = e3.join(e6, on="entryId", how="semi").height
m3 = v3.join(c3, on="repId", how="semi")  # members of v3 non-singleton clusters
res["v3_cluster_members"] = m3.height
res["v3_cluster_members_in_v6"] = m3.select("entryId").join(e6, on="entryId", how="semi").height
res["v3_reps_in_v6_entries"] = c3.join(e6.rename({"entryId": "repId"}), on="repId", how="semi").height
res["v3_reps_also_v6_reps"] = c3.join(c6, on="repId", how="semi").height
# v6 taxId-missing entries: can v3 fill them?
miss = v6.filter(pl.col("taxId") == "0").select("entryId", "repId").unique("entryId")
res["v6_taxId0_entries"] = miss.height
res["v6_taxId0_fillable_from_v3"] = miss.join(v3.filter(pl.col("taxId") != "0").select("entryId"), on="entryId", how="semi").height
dup = v6.group_by("entryId").len().filter(pl.col("len") > 1)
res["v6_entries_listed_more_than_once"] = dup.height
print(json.dumps(res, indent=1))
(ROOT / "results/data_scout").mkdir(parents=True, exist_ok=True)
json.dump(res, open(ROOT / "results/data_scout/afdb_v3_v6_overlap.json", "w"), indent=1)
