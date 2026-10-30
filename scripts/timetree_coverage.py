#!/usr/bin/env python3
"""Coverage of AFDB (v3 clusters) species by the TimeTree of Life 2015 species tree (Hedges et al. 2015).

Maps TimeTree tip names (Genus_species) to NCBI taxIds via scientific names in names.dmp,
then counts how many AFDB proteins / species / clusters fall in species present in the tree.
Output: results/data_scout/timetree_coverage.json, data/interim/timetree/ttol2015_tips_taxid.tsv
"""
import json
import re
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[1]
TT = ROOT / "data/raw/timetree/TimetreeOfLife2015.nwk"
TAX = ROOT / "data/raw/ncbi_taxonomy/new_taxdump"
OUTI = ROOT / "data/interim/timetree"
OUTI.mkdir(parents=True, exist_ok=True)
OUT = ROOT / "results/data_scout"


def dmp(path, n):
    df = pl.read_csv(path, separator="\t", has_header=False, quote_char=None, infer_schema=False)
    df = df.select([df.columns[i] for i in range(0, len(df.columns), 2)][:n])
    return df


s = TT.read_text()
tips = re.findall(r"[(,]([^():,;]+):", s)
tipdf = pl.DataFrame({"tip": tips}).with_columns(pl.col("tip").str.replace_all("_", " ").alias("name"))

names = dmp(TAX / "names.dmp", 4)
names.columns = ["taxId", "name", "uniq", "cls"]
sci = names.filter(pl.col("cls") == "scientific name").select(pl.col("taxId").cast(pl.Int64), "name")
syn = names.filter(pl.col("cls").is_in(["synonym", "equivalent name", "genbank synonym"])).select(
    pl.col("taxId").cast(pl.Int64).alias("taxId_syn"), "name").unique("name", keep="first")
rl = dmp(TAX / "rankedlineage.dmp", 10)
rl.columns = ["taxId", "tname", "species", "genus", "family", "order", "class", "phylum", "kingdom", "domain"]
rl = rl.with_columns(pl.col("taxId").cast(pl.Int64))

m = tipdf.join(sci, on="name", how="left").join(syn, on="name", how="left").with_columns(
    pl.coalesce("taxId", "taxId_syn").alias("taxId")).drop("taxId_syn")
m = m.join(rl.select("taxId", "domain", "kingdom"), on="taxId", how="left")
m.write_csv(OUTI / "ttol2015_tips_taxid.tsv", separator="\t")
res = {"ttol2015_tips": len(tips), "ttol2015_tips_mapped_to_ncbi": int(m["taxId"].is_not_null().sum())}
res["ttol2015_tips_by_domain"] = {str(k): v for k, v in m.group_by("domain").len().sort("len", descending=True).iter_rows()}
res["ttol2015_tips_by_kingdom"] = {str(k): v for k, v in m.group_by("kingdom").len().sort("len", descending=True).head(15).iter_rows()}

# AFDB v3 members: map taxId -> species-level taxId (taxId itself if species; else species name lookup)
v3 = pl.read_csv(ROOT / "data/raw/afdb_clusters/v3/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz",
                 separator="\t", has_header=False, new_columns=["repId", "entryId", "cluFlag", "taxId"],
                 infer_schema_length=0).with_columns(pl.col("taxId").cast(pl.Int64, strict=False))
cnt = v3.group_by("taxId").len()
cnt = cnt.join(rl.select("taxId", "tname", "species", "domain"), on="taxId", how="left").with_columns(
    pl.when(pl.col("species").is_null() | (pl.col("species") == "")).then(pl.col("tname")).otherwise(pl.col("species")).alias("species_name"))
tt_names = set(m.filter(pl.col("taxId").is_not_null())["name"].to_list()) | set(tipdf["name"].to_list())
cnt = cnt.with_columns(pl.col("species_name").is_in(list(tt_names)).alias("in_ttol"))
sp = cnt.group_by("species_name", "domain").agg(pl.col("len").sum(), pl.col("in_ttol").first())
res["afdb_v3_species"] = sp.height
res["afdb_v3_species_in_ttol2015"] = int(sp["in_ttol"].sum())
res["afdb_v3_proteins_in_ttol2015_species"] = int(cnt.filter(pl.col("in_ttol"))["len"].sum())
res["afdb_v3_proteins_total"] = int(cnt["len"].sum())
res["afdb_v3_species_in_ttol_by_domain"] = {str(d): [int(a), int(b)] for d, a, b in sp.group_by("domain").agg(pl.col("in_ttol").sum().alias("a"), pl.len().alias("b")).sort("b", descending=True).iter_rows()}
# clusters with >=2 TTOL species
v3t = v3.join(cnt.select("taxId", "in_ttol", "species_name"), on="taxId", how="left").filter(pl.col("in_ttol"))
c2 = pl.read_csv(ROOT / "data/raw/afdb_clusters/v3/2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv.gz",
                 separator="\t", has_header=False, infer_schema_length=0).select(pl.col("column_0").alias("repId"))
v3t = v3t.join(c2, on="repId", how="semi")
per = v3t.group_by("repId").agg(pl.col("species_name").n_unique().alias("n_tt_species"))
res["v3_clusters_with_ge2_ttol_species"] = int((per["n_tt_species"] >= 2).sum())
res["v3_clusters_with_ge10_ttol_species"] = int((per["n_tt_species"] >= 10).sum())
json.dump(res, open(OUT / "timetree_coverage.json", "w"), indent=1)
print(json.dumps(res, indent=1))
