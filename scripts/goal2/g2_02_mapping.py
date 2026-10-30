#!/usr/bin/env python3
"""Goal 2-ii items 1 and 3: HOG <-> AFDB v3 structural cluster mapping, and cluster-LCA vs HOG-root ages.

Inputs: data/interim/goal2/hog_members.tsv (g2_01), protein_annot.parquet, Ancestra guide tree, NCBI taxdump,
panel-tree node ages (data/interim/timetree; ages stay in data/interim/goal2, never in results/).

Eligible member (primary): AFDB v4 model with md5-identical sequence (seq_match), not a fragment (cluFlag 3),
in a non-singleton AFDB v3 Foldseek cluster (afdb_in_cluster). Cluster id = afdb_cluster_rep.
Sensitivity 'withsingletons': non-clustered, non-fragment entries keep their own structural unit
(afdb_rep_raw = AFDB50 rep that is a Foldseek singleton).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from g2_common import (EUK, OUT_F, OUT_I, OUT_R, ROOT, PanelAges, guide_tree, leaf_taxid, lineage,
                       load_taxonomy)

OUT_R.mkdir(parents=True, exist_ok=True); OUT_F.mkdir(parents=True, exist_ok=True)
summary = {}

mem = pl.read_csv(OUT_I / "hog_members.tsv", separator="\t", schema_overrides={"hog_level": pl.Utf8})
ann = pl.read_parquet(ROOT / "data/processed/annotations/protein_annot.parquet",
                      columns=["acc", "leaf", "panel", "seq_match", "afdb_cluflag", "afdb_in_cluster",
                               "afdb_cluster_rep", "afdb_rep_raw", "afdb_cluster_nmem", "afdb_lca_taxid",
                               "afdb_lca_name", "afdb_lca_group", "afdb_dark"]).filter(pl.col("panel") == "euk")
summary["n_top_hogs"] = mem["hog"].n_unique()
summary["n_genes_in_hogs"] = mem.height
summary["n_euk_panel_proteins"] = ann.height
assert summary["n_top_hogs"] == 109598, summary
m = mem.join(ann, on="acc", how="left")
assert m["panel"].null_count() == 0, "HOG member missing from annotation table"
assert (m["species"] == m["leaf"]).all()

elig = (pl.col("seq_match") & (pl.col("afdb_cluflag") != 3))
m = m.with_columns(
    eligible=(elig & pl.col("afdb_in_cluster")).fill_null(False),
    clu=pl.when(elig & pl.col("afdb_in_cluster")).then(pl.col("afdb_cluster_rep")),
    clu_s=pl.when(elig & pl.col("afdb_cluflag").is_not_null())
            .then(pl.coalesce("afdb_cluster_rep", "afdb_rep_raw")),
)
summary["member_status"] = {
    "eligible_in_cluster": int(m["eligible"].sum()),
    "no_model_or_seq_mismatch": int((~m["seq_match"].fill_null(False)).sum()),
    "fragment": int((m["seq_match"].fill_null(False) & (m["afdb_cluflag"] == 3)).sum()),
    "structural_singleton": int((m["seq_match"].fill_null(False) & (m["afdb_cluflag"] != 3)
                                 & ~m["afdb_in_cluster"]).sum()),
}
m.select("acc", "species", "hog", "hog_level", "hog_size", "eligible", "clu", "clu_s").write_parquet(
    OUT_I / "hog_members_clu.parquet")


def contingency(m, col, tag):
    e = m.filter(pl.col(col).is_not_null())
    # per HOG
    ph = (e.group_by("hog").agg(n_elig=pl.len(), n_clu=pl.col(col).n_unique(),
                                 top=pl.col(col).value_counts(sort=True).first().struct.field("count"))
          .with_columns(frac_largest=pl.col("top") / pl.col("n_elig")))
    # per cluster (panel members that are in HOGs)
    pc = e.group_by(col).agg(n_mem=pl.len(), n_hog=pl.col("hog").n_unique())
    # flags
    hog_clus = e.select("hog", col).unique()
    multi_clu = pc.filter(pl.col("n_hog") > 1).select(col)
    merged_h = hog_clus.join(multi_clu, on=col).select("hog").unique().with_columns(merged=pl.lit(True))
    ph = ph.join(merged_h, on="hog", how="left").with_columns(pl.col("merged").fill_null(False))
    ph2 = ph.filter(pl.col("n_elig") >= 2).with_columns(split=pl.col("n_clu") > 1)
    cat = (pl.when(~pl.col("split") & ~pl.col("merged")).then(pl.lit("1HOG_1cluster"))
           .when(pl.col("split") & ~pl.col("merged")).then(pl.lit("HOG_split_only"))
           .when(~pl.col("split") & pl.col("merged")).then(pl.lit("cluster_merges_only"))
           .otherwise(pl.lit("split_and_merged")))
    ph2 = ph2.with_columns(category=cat)
    tab = ph2.group_by("category").agg(n_hogs=pl.len(), n_members=pl.col("n_elig").sum()).sort("category")
    tab = tab.with_columns(frac_hogs=pl.col("n_hogs") / pl.col("n_hogs").sum(),
                           frac_members=pl.col("n_members") / pl.col("n_members").sum())
    # cluster-level categories (clusters with >= 2 HOG-assigned panel members)
    hog_nclu = ph.select("hog", "n_clu")
    cl = (e.join(hog_nclu, on="hog").group_by(col)
          .agg(n_mem=pl.len(), n_hog=pl.col("hog").n_unique(), any_split=(pl.col("n_clu") > 1).any()))
    cl2 = cl.filter(pl.col("n_mem") >= 2)
    ccat = (pl.when((pl.col("n_hog") == 1) & ~pl.col("any_split")).then(pl.lit("1cluster_1HOG"))
            .when((pl.col("n_hog") == 1) & pl.col("any_split")).then(pl.lit("cluster_is_part_of_split_HOG"))
            .when((pl.col("n_hog") > 1) & ~pl.col("any_split")).then(pl.lit("cluster_merges_whole_HOGs"))
            .otherwise(pl.lit("merges_and_part_of_split")))
    ctab = (cl2.with_columns(category=ccat).group_by("category").agg(n_clusters=pl.len(), n_members=pl.col("n_mem").sum())
            .sort("category").with_columns(frac_clusters=pl.col("n_clusters") / pl.col("n_clusters").sum()))
    tab.write_csv(OUT_R / f"contingency_hog_{tag}.tsv", separator="\t")
    ctab.write_csv(OUT_R / f"contingency_cluster_{tag}.tsv", separator="\t")
    s = {
        "n_hogs_with_ge1_eligible": ph.height,
        "n_hogs_with_ge2_eligible": ph2.height,
        "hog_n_clusters_quantiles": {q: float(ph2["n_clu"].quantile(q)) for q in (0.5, 0.9, 0.99)},
        "hog_frac_split": float(ph2["split"].mean()),
        "hog_frac_largest_median_split": float(ph2.filter(pl.col("split"))["frac_largest"].median()),
        "hog_frac_largest_mean_all": float(ph2["frac_largest"].mean()),
        "member_weighted_frac_in_largest": float(ph2["top"].sum() / ph2["n_elig"].sum()),
        "n_clusters_with_panel_hog_members": cl.height,
        "cluster_n_hogs_quantiles": {q: float(cl["n_hog"].quantile(q)) for q in (0.5, 0.9, 0.99)},
        "cluster_frac_multi_hog_ge2mem": float((cl2["n_hog"] > 1).mean()),
        "hog_table": tab.to_dicts(), "cluster_table": ctab.to_dicts(),
    }
    return ph2, cl, s


ph, cl, summary["primary"] = contingency(m, "clu", "primary")
_, _, summary["sens_with_singletons"] = contingency(m, "clu_s", "withsingletons")
ph.write_parquet(OUT_I / "hog_cluster_stats.parquet"); cl.write_parquet(OUT_I / "cluster_hog_stats.parquet")

# HOG size dependence (size = eligible members) of the split rate
sz = (ph.with_columns(bin=pl.col("n_elig").cut([2, 5, 10, 20, 50, 100, 200, 500],
                                                 labels=["2", "3-5", "6-10", "11-20", "21-50", "51-100", "101-200",
                                                         "201-500", ">500"]))
      .group_by("bin").agg(n=pl.len(), frac_split=pl.col("split").mean(),
                           mean_nclu=pl.col("n_clu").mean(), frac_largest=pl.col("frac_largest").mean()).sort("bin"))
sz.write_csv(OUT_R / "split_by_hog_size.tsv", separator="\t")

# ---------------------------------------------------------------- item 3: ages
par, rank, names = load_taxonomy()
leafsets, leaves, _ = guide_tree()
leaf_lin = {l: lineage(leaf_taxid(l), par) for l in leaves}
leaf_linset = {l: set(v) for l, v in leaf_lin.items()}
pa = PanelAges()
lvl_ok = set(m["hog_level"].unique().to_list()) - set(leafsets)
summary["hog_levels_not_in_guide_tree"] = sorted(map(str, lvl_ok))[:20]
root_name = "Eukaryota"
root_leafset = leafsets[root_name]
assert len(root_leafset) == 180


def proj_leafset(t):
    """Panel leaves descending from NCBI taxon t."""
    return frozenset(l for l in leaves if t in leaf_linset[l])


def is_euk(t):
    return EUK in set(lineage(t, par))


# per HOG: root-level leaf set; panel-restricted HOG: leaf set of species present
hog_info = (m.group_by("hog").agg(level=pl.col("hog_level").first(), kind=pl.col("hog_kind").first(),
                                  species=pl.col("species").unique()))
hog_lv = {}
for r in hog_info.iter_rows(named=True):
    lv = leafsets.get(str(r["level"]))
    if lv is None:  # leaf-level HOG (species-specific) or unknown name
        lv = frozenset([str(r["level"])]) if str(r["level"]) in leaf_lin else frozenset(r["species"])
    hog_lv[r["hog"]] = lv

# per cluster with eligible HOG members: full-AFDB LCA, panel-restricted LCA (all euk panel members of the
# cluster, in a HOG or not), dominant HOG (most members; ties -> oldest = largest leafset)
e = m.filter(pl.col("clu").is_not_null())
cl_all = (ann.filter(pl.col("afdb_in_cluster") & pl.col("seq_match") & (pl.col("afdb_cluflag") != 3))
          .group_by("afdb_cluster_rep").agg(lca=pl.col("afdb_lca_taxid").first(), lca_group=pl.col("afdb_lca_group").first(),
                                            nmem=pl.col("afdb_cluster_nmem").first(), panel_leaves=pl.col("leaf").unique()))
dom = (e.group_by("clu", "hog").agg(n=pl.len()).sort(["clu", "n"], descending=[False, True])
       .group_by("clu").agg(dom_hog=pl.col("hog").first(), dom_n=pl.col("n").first(), n_hog=pl.len(),
                             hogs=pl.col("hog")))
cdf = dom.join(cl_all.rename({"afdb_cluster_rep": "clu"}), on="clu", how="left")

lca_cache = {}


def cmp_sets(a, b):
    """a, b: leaf sets on the same lineage (nested). Returns equal / a_older / a_younger / incomparable."""
    if a == b:
        return "equal"
    if b < a:
        return "a_older"
    if a < b:
        return "a_younger"
    return "incomparable"


rows = []
for r in cdf.iter_rows(named=True):
    t = r["lca"]
    if t is None:
        continue
    t = int(t)
    if t not in lca_cache:
        if is_euk(t):
            lca_cache[t] = proj_leafset(t)
        else:
            lca_cache[t] = None
    cl_proj = lca_cache[t]
    # panel-restricted cluster LCA = guide-tree MRCA of panel species (as leaf set: minimal guide node containing)
    pl_leaves = frozenset(r["panel_leaves"])
    # smallest guide-tree clade containing pl_leaves
    best = None
    for nm, lv in leafsets.items():
        if pl_leaves <= lv and (best is None or len(lv) < len(best)):
            best = lv
    hog_set = hog_lv[r["dom_hog"]]
    # oldest HOG among all HOGs in the cluster
    oldest = max((hog_lv[h] for h in r["hogs"]), key=len)
    full_vs_hog = "cluster_pre_euk" if cl_proj is None else cmp_sets(cl_proj, hog_set)
    rows.append(dict(
        clu=r["clu"], n_hog=r["n_hog"], dom_hog=r["dom_hog"], nmem=r["nmem"], lca_taxid=t,
        lca_name=names.get(t, ""), lca_group=r["lca_group"],
        lca_pre_euk=cl_proj is None,
        hog_at_root=hog_set == root_leafset,
        full_vs_domhog=full_vs_hog,
        full_vs_oldesthog="cluster_pre_euk" if cl_proj is None else cmp_sets(cl_proj, oldest),
        panel_vs_domhog=cmp_sets(best, hog_set),
        full_vs_panel="pre_euk" if cl_proj is None else cmp_sets(cl_proj, best),
        hog_age=pa.mrca(hog_set)[1], panel_lca_age=pa.mrca(best)[1],
        full_lca_age=(float("inf") if cl_proj is None else pa.mrca(cl_proj)[1]),
    ))
ages = pl.DataFrame(rows)
ages.write_parquet(OUT_I / "cluster_lca_vs_hog_ages.parquet")  # ages: interim only


def consistency(df):
    """Per cluster: consistent / cluster older / cluster younger, with explicit censoring at the euk root."""
    c = (pl.when(pl.col("full_vs_domhog") == "equal").then(pl.lit("consistent"))
         .when((pl.col("full_vs_domhog") == "cluster_pre_euk") & pl.col("hog_at_root")).then(pl.lit("consistent_censored_pre_euk"))
         .when(pl.col("full_vs_domhog") == "cluster_pre_euk").then(pl.lit("cluster_older_pre_euk"))
         .when(pl.col("full_vs_domhog") == "a_older").then(pl.lit("cluster_older_within_euk"))
         .when(pl.col("full_vs_domhog") == "a_younger").then(pl.lit("cluster_younger"))
         .otherwise(pl.lit("incomparable")))
    return df.with_columns(cons=c)


ages = consistency(ages)
ctab = (ages.group_by("cons").agg(n=pl.len()).with_columns(frac=pl.col("n") / pl.col("n").sum()).sort("cons"))
ctab.write_csv(OUT_R / "lca_consistency_clusters.tsv", separator="\t")
# decomposition: full-AFDB LCA vs panel-restricted LCA (sampling) and panel LCA vs HOG (partition)
dec = ages.group_by("full_vs_panel", "panel_vs_domhog").agg(n=pl.len()).with_columns(
    frac=pl.col("n") / pl.col("n").sum()).sort("n", descending=True)
dec.write_csv(OUT_R / "lca_decomposition_clusters.tsv", separator="\t")
# by HOG root-age stratum (labels only; no ages written to results)
summary["lca_consistency"] = ctab.to_dicts()
summary["lca_consistency_euk_lca_only"] = (
    ages.filter(~pl.col("lca_pre_euk")).group_by("cons").agg(n=pl.len())
    .with_columns(frac=pl.col("n") / pl.col("n").sum()).sort("cons").to_dicts())
summary["lca_consistency_by_cluster_nhog"] = (
    ages.with_columns(multi=pl.col("n_hog") > 1).group_by("multi", "cons").agg(n=pl.len())
    .with_columns(frac=pl.col("n") / pl.col("n").sum().over("multi")).sort("multi", "cons").to_dicts())
# per-HOG view: HOG root vs LCA of its largest cluster
hl = (e.group_by("hog", "clu").agg(n=pl.len()).sort(["hog", "n"], descending=[False, True])
      .group_by("hog").agg(clu=pl.col("clu").first()))
hv = hl.join(ages.select("clu", "lca_pre_euk", "lca_taxid"), on="clu", how="left")
hv_rows = []
for r in hv.iter_rows(named=True):
    if r["lca_taxid"] is None:
        continue
    hs = hog_lv[r["hog"]]
    if r["lca_pre_euk"]:
        c = "consistent_censored_pre_euk" if hs == root_leafset else "cluster_older_pre_euk"
    else:
        k = cmp_sets(lca_cache[r["lca_taxid"]], hs)
        c = {"equal": "consistent", "a_older": "cluster_older_within_euk", "a_younger": "cluster_younger"}.get(k, k)
    hv_rows.append((r["hog"], c, hs == root_leafset))
hv = pl.DataFrame(hv_rows, schema=["hog", "cons", "hog_at_root"], orient="row")
htab = hv.group_by("cons").agg(n=pl.len()).with_columns(frac=pl.col("n") / pl.col("n").sum()).sort("cons")
htab.write_csv(OUT_R / "lca_consistency_hogs.tsv", separator="\t")
summary["lca_consistency_hogs"] = htab.to_dicts()
summary["hog_root_at_euk_root_frac"] = float(hv["hog_at_root"].mean())
json.dump(summary, open(OUT_R / "mapping_summary.json", "w"), indent=1, default=str)
print(json.dumps(summary, indent=1, default=str)[:6000])
