#!/usr/bin/env python3
"""Quick descriptive statistics of the AFDB Clusters tables (v3 and v6 releases).

Inputs (downloaded by scripts/download_afdb_clusters.sh, NCBI taxdump in data/raw/ncbi_taxonomy):
  data/raw/afdb_clusters/{v3,v6}/2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv.gz
  data/raw/afdb_clusters/{v3,v6}/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz
  data/raw/ncbi_taxonomy/new_taxdump/{rankedlineage,taxidlineage,merged}.dmp
Outputs: results/data_scout/afdb_cluster_stats_{v3,v6}.json and per-cluster taxonomic span tables (parquet).

Usage: python scripts/afdb_cluster_stats.py v3|v6
"""
import json
import sys
import time
from pathlib import Path

import polars as pl

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/afdb_clusters"
TAX = ROOT / "data/raw/ncbi_taxonomy/new_taxdump"
OUT = ROOT / "results/data_scout"
OUT.mkdir(parents=True, exist_ok=True)

# AFDB model-organism proteome taxIds (as used in AFDB proteome downloads) -> label
MODEL_ORGS = {
    9606: "H.sapiens", 10090: "M.musculus", 10116: "R.norvegicus", 7955: "D.rerio",
    7227: "D.melanogaster", 6239: "C.elegans", 3702: "A.thaliana", 559292: "S.cerevisiae",
    284812: "S.pombe", 44689: "D.discoideum", 83333: "E.coli_K12", 243232: "M.jannaschii",
    237561: "C.albicans", 4577: "Z.mays", 39947: "O.sativa", 3847: "G.max",
}
DOMAINS = {2: "Bacteria", 2157: "Archaea", 2759: "Eukaryota", 10239: "Viruses"}


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def read_dmp(path, cols):
    # NCBI dmp: fields separated by "\t|\t", rows end with "\t|"
    df = pl.read_csv(path, separator="\t", has_header=False, quote_char=None,
                     infer_schema=False)
    # keep every other column (the "|" separators are their own columns)
    df = df.select([df.columns[i] for i in range(0, len(df.columns), 2)][: len(cols)])
    df.columns = cols
    return df


def taxonomy_table():
    """taxId -> domain, kingdom, phylum, species_taxid-free labels; resolve merged ids."""
    log("loading taxonomy")
    rl = read_dmp(TAX / "rankedlineage.dmp",
                  ["taxId", "name", "species", "genus", "family", "order", "class",
                   "phylum", "kingdom", "domain"])
    tl = read_dmp(TAX / "taxidlineage.dmp", ["taxId", "lineage"])
    rl = rl.with_columns(pl.col("taxId").cast(pl.Int64))
    tl = tl.with_columns(pl.col("taxId").cast(pl.Int64))
    # domain via lineage membership (robust for the domain nodes themselves and viruses)
    lin = tl.with_columns(pl.col("lineage").str.strip_chars().str.split(" ").alias("l"))
    expr = None
    for tid, name in DOMAINS.items():
        cond = pl.col("l").list.contains(str(tid)) | (pl.col("taxId") == tid)
        expr = pl.when(cond).then(pl.lit(name)) if expr is None else expr.when(cond).then(pl.lit(name))
    lin = lin.with_columns(expr.otherwise(pl.lit("Other")).alias("domain2")).select("taxId", "domain2")
    tax = rl.join(lin, on="taxId", how="left").select(
        "taxId", "name", "species", "phylum", "kingdom", pl.col("domain2").alias("domain"))
    tax = tax.with_columns(
        pl.when(pl.col("kingdom").is_null() | (pl.col("kingdom") == ""))
        .then(pl.col("domain")).otherwise(pl.col("kingdom")).alias("kingdom_or_domain"))
    merged = read_dmp(TAX / "merged.dmp", ["old", "new"]).with_columns(
        pl.col("old").cast(pl.Int64), pl.col("new").cast(pl.Int64))
    return tax, merged


def resolve(df, tax, merged, col="taxId"):
    df = df.join(merged, left_on=col, right_on="old", how="left").with_columns(
        pl.coalesce(pl.col("new"), pl.col(col)).alias(col + "_res")).drop("new")
    return df.join(tax.select("taxId", "domain", "kingdom_or_domain", "phylum"),
                   left_on=col + "_res", right_on="taxId", how="left")


def main(ver):
    res = {"version": ver}
    tax, merged = taxonomy_table()

    # ---------------- file 2: cluster overview ----------------
    f2 = RAW / ver / "2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv.gz"
    c = pl.read_csv(f2, separator="\t", has_header=False, infer_schema_length=0,
                    new_columns=["repId", "isDark", "nMem", "repLen", "avgLen", "repPlddt",
                                 "avgPlddt", "LCAtaxId"])
    c = c.with_columns(pl.col("isDark").cast(pl.Int8), pl.col("nMem").cast(pl.Int64),
                       pl.col("repLen").cast(pl.Int64), pl.col("avgLen").cast(pl.Float64),
                       pl.col("repPlddt").cast(pl.Float64), pl.col("avgPlddt").cast(pl.Float64),
                       pl.col("LCAtaxId").cast(pl.Int64, strict=False))
    log("file2 rows", c.height)
    res["file2_n_clusters"] = c.height
    res["file2_n_dark"] = int(c["isDark"].sum())
    res["file2_sum_nMem"] = int(c["nMem"].sum())
    q = [0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 0.999, 1.0]
    res["file2_nMem_quantiles"] = {str(x): float(c["nMem"].quantile(x)) for x in q}
    bins = [1, 2, 5, 10, 20, 50, 100, 1000, 10000, 10**9]
    res["file2_nMem_bins"] = {
        f"[{bins[i]},{bins[i+1]})": int(c.filter((pl.col("nMem") >= bins[i]) & (pl.col("nMem") < bins[i + 1])).height)
        for i in range(len(bins) - 1)}
    res["file2_avgPlddt_quantiles"] = {str(x): float(c["avgPlddt"].quantile(x)) for x in q}
    res["file2_n_avgPlddt_ge70"] = int((c["avgPlddt"] >= 70).sum())
    res["file2_n_avgPlddt_ge90"] = int((c["avgPlddt"] >= 90).sum())
    cl = resolve(c, tax, merged, "LCAtaxId")
    res["file2_LCA_domain_counts"] = {str(k): v for k, v in cl.group_by("domain").len().sort("len", descending=True).iter_rows()}
    lca_names = cl.join(tax.select("taxId", "name"), left_on="LCAtaxId_res", right_on="taxId", how="left")
    res["file2_LCA_top_taxa"] = [[n, t, int(k)] for n, t, k in lca_names.group_by("name", "LCAtaxId").len().sort("len", descending=True).head(25).iter_rows()]

    # ---------------- file 5: all members ----------------
    f5 = RAW / ver / "5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz"
    log("reading file5")
    m = pl.read_csv(f5, separator="\t", has_header=False, infer_schema_length=0,
                    new_columns=["repId", "entryId", "cluFlag", "taxId"])
    m = m.with_columns(pl.col("cluFlag").cast(pl.Int8), pl.col("taxId").cast(pl.Int64, strict=False))
    log("file5 rows", m.height)
    res["file5_n_rows"] = m.height
    res["file5_n_unique_entries"] = int(m["entryId"].n_unique())
    res["file5_n_unique_reps"] = int(m["repId"].n_unique())
    res["file5_cluFlag_counts"] = {str(k): v for k, v in m.group_by("cluFlag").len().sort("cluFlag").iter_rows()}
    res["file5_n_taxId0_or_null"] = int(((m["taxId"] == 0) | m["taxId"].is_null()).sum())
    res["file5_n_unique_taxIds"] = int(m["taxId"].n_unique())

    # members belonging to the non-singleton, fragment-free AFDB clusters (reps in file 2)
    reps = c.select("repId")
    mc = m.join(reps, on="repId", how="semi")
    log("members of file2 clusters", mc.height)
    res["members_in_file2_clusters_n"] = mc.height
    res["members_in_file2_clusters_cluFlag"] = {str(k): v for k, v in mc.group_by("cluFlag").len().sort("cluFlag").iter_rows()}
    mc = resolve(mc, tax, merged, "taxId")
    res["members_domain_counts"] = {str(k): v for k, v in mc.group_by("domain").len().sort("len", descending=True).iter_rows()}

    # per-cluster taxonomic span (all members incl. AFDB50 members)
    span = mc.group_by("repId").agg(
        pl.len().alias("n_all_members"),
        pl.col("taxId_res").n_unique().alias("n_taxa"),
        pl.col("domain").filter(pl.col("domain").is_in(["Bacteria", "Archaea", "Eukaryota"])).n_unique().alias("n_cell_domains"),
        pl.col("domain").filter(pl.col("domain").is_in(["Bacteria", "Archaea", "Eukaryota", "Viruses"])).n_unique().alias("n_domains_incl_virus"),
        pl.col("kingdom_or_domain").filter(pl.col("kingdom_or_domain").is_not_null() & (pl.col("kingdom_or_domain") != "Other")).n_unique().alias("n_kingdoms"),
        pl.col("phylum").filter(pl.col("phylum").is_not_null() & (pl.col("phylum") != "")).n_unique().alias("n_phyla"),
        pl.col("domain").filter(pl.col("domain").is_not_null()).unique().sort().str.join(";").alias("domains"),
        *[(pl.col("taxId_res") == t).sum().alias(f"n_{lab}") for t, lab in MODEL_ORGS.items()],
    )
    span = span.join(c, on="repId", how="left")
    span.write_parquet(OUT / f"afdb_clusters_{ver}_taxspan.parquet")
    res["span_n_clusters"] = span.height
    res["span_n_clusters_ge2_cell_domains"] = int((span["n_cell_domains"] >= 2).sum())
    res["span_n_clusters_all3_cell_domains"] = int((span["n_cell_domains"] == 3).sum())
    res["span_n_clusters_ge2_kingdoms"] = int((span["n_kingdoms"] >= 2).sum())
    res["span_n_clusters_ge5_phyla"] = int((span["n_phyla"] >= 5).sum())
    res["span_domain_combos"] = {str(k): v for k, v in span.group_by("domains").len().sort("len", descending=True).iter_rows()}
    res["span_n_taxa_quantiles"] = {str(x): float(span["n_taxa"].quantile(x)) for x in q}
    mo = {}
    for t, lab in MODEL_ORGS.items():
        col = f"n_{lab}"
        mo[lab] = {"taxId": t, "clusters_with_member": int((span[col] > 0).sum()),
                   "proteins_in_clusters": int(span[col].sum()),
                   "proteins_in_file5": int((m["taxId"] == t).sum())}
    res["model_organisms"] = mo
    mcols = [f"n_{lab}" for lab in MODEL_ORGS.values()]
    res["span_n_clusters_any_model_org"] = int(span.select(pl.any_horizontal([pl.col(x) > 0 for x in mcols])).to_series().sum())
    res["span_n_clusters_human_and_ecoli"] = int(((span["n_H.sapiens"] > 0) & (span["n_E.coli_K12"] > 0)).sum())
    res["span_n_clusters_human_and_yeast"] = int(((span["n_H.sapiens"] > 0) & (span["n_S.cerevisiae"] > 0)).sum())
    res["span_n_clusters_ge3_model_orgs"] = int(span.select(pl.sum_horizontal([(pl.col(x) > 0).cast(pl.Int32) for x in mcols]) >= 3).to_series().sum())

    with open(OUT / f"afdb_cluster_stats_{ver}.json", "w") as fh:
        json.dump(res, fh, indent=1)
    log("done")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
