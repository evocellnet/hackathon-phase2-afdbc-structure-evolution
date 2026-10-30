#!/usr/bin/env python3
"""Assemble per-protein annotation for all panel proteins -> data/processed/annotations/protein_annot.parquet
and results/annotations/summary.json.  Inputs: data/interim/annotations/{afdb_cluster_annot,ted_domains_panel,ted_chains_panel,pfam_regions_panel}.parquet,
panel protein tables, pfamA.txt.gz (family type).  Run after afdb_cluster_annot.py, ted_extract.py {domains,chains}, pfam_extract.py.

Class rule (ANALYSIS_PLAN 4.3): CATH class digit (1 alpha, 2 beta, 3 alpha/beta) covering >= 60% of TED-domain residues (all TED domains of the protein,
high+medium consensus, domains lacking a CATH label stay in the denominator); else `mixed`; `none` if the v4 model has no TED domain; `no_model` if no AFDB v4 model.
If class 4 (few secondary structures) covers >= 60% the protein is `mixed` (not alpha/beta/alphabeta).  `ted_class_labeled` is a sensitivity version whose
denominator is only CATH-labelled residues (`unlabeled` if no domain has a CATH label). Note: under the literal rule, proteins whose TED domains all lack a CATH label are `mixed` in `ted_class`.
"""
import json, gzip
from pathlib import Path
import polars as pl
ROOT = Path(__file__).resolve().parents[2]
I = ROOT/"data/interim/annotations"; O = ROOT/"data/processed/annotations"; O.mkdir(parents=True, exist_ok=True)
R = ROOT/"results/annotations"; R.mkdir(parents=True, exist_ok=True)

pt = pl.concat([pl.read_parquet(ROOT/f"data/processed/panel/{f}.parquet").with_columns(panel=pl.lit(p))
                for f, p in [("protein_table", "euk"), ("protein_table_prok", "prok")]])
base = pt.select("acc", "upid", "leaf", "name", "superkingdom", "panel", "has_v4", "seq_match", "seq_match_ted", "len_2026", "len_v4", "md5_2026", "mean_plddt")

# ---------- 1. AFDB clusters ----------
clu = pl.read_parquet(I/"afdb_cluster_annot.parquet")

# ---------- 2. TED ----------
d = pl.read_parquet(I/"ted_domains_panel.parquet")
d = d.with_columns(
    pl.col("chopping").str.extract(r"^(\d+)", 1).cast(pl.Int32).alias("start"),
    pl.col("chopping").str.extract(r"-(\d+)$", 1).cast(pl.Int32).alias("end"),
    pl.when(pl.col("cath_label").is_in(["-", ""]) | pl.col("cath_label").is_null()).then(None).otherwise(pl.col("cath_label")).alias("cath"),
).with_columns(
    pl.col("cath").str.split(".").list.first().alias("cdig"),
    pl.col("cath").str.count_matches(r"\.").add(1).alias("clevels"),  # 4 = C.A.T.H, 3 = C.A.T
).sort("acc", "start")
dom = d.group_by("acc", maintain_order=True).agg(
    pl.len().alias("ted_n_domains"),
    (pl.col("consensus") == "high").sum().alias("ted_n_high"),
    (pl.col("consensus") == "medium").sum().alias("ted_n_medium"),
    pl.col("cath").is_not_null().sum().alias("ted_n_cath"),
    (pl.col("clevels") == 4).sum().alias("ted_n_cath_H"),
    pl.col("chopping").str.join(";").alias("ted_boundaries"),
    pl.col("cath").drop_nulls().str.join("|").alias("ted_arch_cath"),
    pl.col("cath").fill_null("NA").str.join("|").alias("ted_arch_all"),
    pl.col("consensus").str.join("|").alias("ted_consensus"),
    pl.col("dom_nres").sum().alias("ted_dom_residues"),
    ((pl.col("dom_plddt") * pl.col("dom_nres")).sum() / pl.col("dom_nres").sum()).alias("ted_mean_dom_plddt"),
    pl.col("dom_plddt").min().alias("ted_min_dom_plddt"),
    *[pl.col("dom_nres").filter(pl.col("cdig") == k).sum().alias(f"_r{k}") for k in "1234"],
    pl.col("dom_nres").filter(pl.col("cdig").is_not_null() & ~pl.col("cdig").is_in(list("1234"))).sum().alias("_r_other"),
    pl.col("dom_nres").filter(pl.col("cdig").is_null()).sum().alias("_r_unl"),
    pl.col("cdig").drop_nulls().unique().sort().str.join("").alias("ted_class_digits_present"),
)
dom = dom.with_columns(pl.col("ted_arch_cath").replace("", None))
dom = dom.with_columns(
    (pl.col("_r1") + pl.col("_r2") + pl.col("_r3") + pl.col("_r4") + pl.col("_r_other")).alias("_rlab"),
)
for k in "1234":
    dom = dom.with_columns((pl.col(f"_r{k}") / pl.col("ted_dom_residues")).alias(f"ted_frac_class{k}"))
dom = dom.with_columns(
    pl.max_horizontal("_r1", "_r2", "_r3", "_r4", "_r_other").alias("_rmax"),
    (pl.col("_r_unl") / pl.col("ted_dom_residues")).alias("ted_frac_unlabeled"),
)
def cls(frac_expr_den):
    out = pl.lit("mixed")
    for k, nm in [("3", "alphabeta"), ("2", "beta"), ("1", "alpha")]:
        out = pl.when(pl.col(f"_r{k}") >= 0.6 * frac_expr_den).then(pl.lit(nm)).otherwise(out)
    return out
dom = dom.with_columns(
    cls(pl.col("ted_dom_residues")).alias("ted_class"),
    pl.when(pl.col("_rlab") > 0).then(cls(pl.col("_rlab"))).otherwise(pl.lit("unlabeled")).alias("ted_class_labeled"),
    pl.when(pl.col("_rlab") > 0).then(
        pl.concat_str([pl.when(pl.col("_r1") == pl.col("_rmax")).then(pl.lit("1")).when(pl.col("_r2") == pl.col("_rmax")).then(pl.lit("2"))
                       .when(pl.col("_r3") == pl.col("_rmax")).then(pl.lit("3")).when(pl.col("_r4") == pl.col("_rmax")).then(pl.lit("4")).otherwise(pl.lit("other"))])
    ).alias("ted_top_class_digit"),
    (pl.col("_rmax") / pl.col("ted_dom_residues")).alias("ted_top_class_frac"),
).drop([c for c in dom.columns if c.startswith("_")])

ch = pl.read_parquet(I/"ted_chains_panel.parquet").select("acc", pl.col("nres").alias("ted_chain_nres"), pl.col("n_low").alias("ted_n_low_consensus_chain"))

# ---------- 3. Pfam ----------
pf = pl.read_parquet(I/"pfam_regions_panel.parquet")
with gzip.open(ROOT/"data/raw/pfam/pfamA.txt.gz", "rt") as fh:
    ptype = pl.DataFrame([(r[0], r[6]) for r in (l.split("\t") for l in fh)], schema=["pfam_acc", "pfam_type"], orient="row")
pf = pf.join(ptype, on="pfam_acc", how="left").join(base.select("acc", "md5_2026"), on="acc", how="left")
pf = pf.with_columns((pl.col("pfam_seq_md5") == pl.col("md5_2026")).alias("md5_ok")).sort("acc", "env_start", "env_end")
# merged envelope coverage
cov = (pf.filter(pl.col("md5_ok")).sort("acc", "env_start").with_columns(
    pl.col("env_end").cum_max().shift(1).over("acc").fill_null(0).alias("_prev_end"))
    .with_columns((pl.col("env_end") - pl.max_horizontal("env_start", pl.col("_prev_end") + 1) + 1).clip(lower_bound=0).alias("_new"))
    .group_by("acc").agg(pl.col("_new").sum().alias("pfam_residues")))
pfa = pf.group_by("acc", maintain_order=True).agg(
    pl.len().alias("pfam_n_regions"),
    pl.col("pfam_acc").str.join("|").alias("pfam_arch"),
    pl.col("pfam_acc").filter(pl.col("pfam_type").is_in(["Domain", "Repeat"])).str.join("|").alias("pfam_arch_domain_repeat"),
    pl.col("pfam_acc").n_unique().alias("pfam_n_families"),
    pl.format("{}-{}", pl.col("env_start"), pl.col("env_end")).str.join(";").alias("pfam_boundaries"),
    pl.col("md5_ok").all().alias("pfam_seq_match_2026"),
).join(cov, on="acc", how="left").with_columns(pl.col("pfam_arch_domain_repeat").replace("", None))

# ---------- assemble ----------
a = (base.join(clu, on="acc", how="left").join(dom, on="acc", how="left").join(ch, on="acc", how="left").join(pfa, on="acc", how="left")
     .join(pf.group_by("acc").agg(pl.lit(True).alias("pfam_has_hit")), on="acc", how="left"))
a = a.with_columns(
    pl.col("pfam_has_hit").fill_null(False),
    pl.col("pfam_n_regions").fill_null(0),
    pl.col("ted_n_domains").fill_null(0),
    (pl.col("ted_n_domains").is_not_null() & (pl.col("ted_n_domains") > 0)).alias("ted_has_domain"),
    pl.col("afdb_in_cluster").fill_null(False),
)
a = a.with_columns(
    pl.when(~pl.col("has_v4")).then(pl.lit("no_model")).when(pl.col("ted_class").is_null()).then(pl.lit("none")).otherwise(pl.col("ted_class")).alias("ted_class"),
    pl.when(~pl.col("has_v4")).then(pl.lit("no_model")).when(pl.col("ted_class_labeled").is_null()).then(pl.lit("none")).otherwise(pl.col("ted_class_labeled")).alias("ted_class_labeled"),
    pl.when(pl.col("ted_dom_residues").is_not_null()).then((pl.col("ted_dom_residues") / pl.col("len_v4")).clip(upper_bound=1.0))
      .when(pl.col("has_v4")).then(pl.lit(0.0)).alias("ted_frac_in_domain"),
    pl.when(pl.col("pfam_residues").is_not_null()).then((pl.col("pfam_residues") / pl.col("len_2026")).clip(upper_bound=1.0)).otherwise(pl.lit(0.0)).alias("pfam_frac_covered"),
    pl.when(~pl.col("has_v4")).then(None).otherwise(pl.col("ted_n_domains")).alias("ted_n_domains"),
)
a = a.drop("md5_2026")
a.write_parquet(O/"protein_annot.parquet", compression="zstd")
print(a.shape); print(a.schema)

# ---------- summary ----------
def stats(x: pl.DataFrame):
    n = x.height; m = x.filter(pl.col("has_v4"))
    nm = m.height
    s = {"n_proteins": n, "n_with_v4_model": nm, "pct_with_v4_model": 100*nm/n}
    s["pct_of_modelled_with_ted_domain"] = 100*m["ted_has_domain"].mean() if nm else None
    s["pct_of_modelled_with_cath_label"] = 100*(m["ted_n_cath"].fill_null(0) > 0).mean() if nm else None
    s["ted_n_domains_mean_modelled"] = m["ted_n_domains"].mean() if nm else None
    s["ted_frac_in_domain_mean_modelled"] = m["ted_frac_in_domain"].mean() if nm else None
    s["class_distribution_pct_of_modelled"] = {k: 100*v/nm for k, v in m.group_by("ted_class").len().iter_rows()} if nm else {}
    ml = m.filter(pl.col("ted_has_domain"))
    s["class_labeled_distribution_pct_of_modelled"] = {k: 100*v/nm for k, v in m.group_by("ted_class_labeled").len().iter_rows()} if nm else {}
    s["pct_of_modelled_ted_proteins_with_all_domains_unlabeled"] = 100*(m["ted_has_domain"] & (m["ted_n_cath"].fill_null(0) == 0)).sum()/nm if nm else None
    s["class_distribution_pct_of_ted_proteins"] = {k: 100*v/ml.height for k, v in ml.group_by("ted_class").len().iter_rows()} if ml.height else {}
    s["pct_of_all_in_afdb_cluster_entry"] = 100*x["afdb_cluflag"].is_not_null().mean()
    inc = x.filter(pl.col("afdb_cluflag").is_not_null())
    s["pct_of_afdb_entries_in_nonsingleton_cluster"] = 100*inc["afdb_in_cluster"].mean()
    s["pct_of_all_in_nonsingleton_cluster"] = 100*x["afdb_in_cluster"].mean()
    s["afdb_cluflag_counts"] = {str(k): v for k, v in inc.group_by("afdb_cluflag").len().sort("afdb_cluflag").iter_rows()}
    c = inc.filter(pl.col("afdb_in_cluster"))
    if c.height:
        s["cluster_nmem_median"] = c["afdb_cluster_nmem"].median(); s["cluster_nmem_p90"] = c["afdb_cluster_nmem"].quantile(0.9)
        s["pct_in_dark_cluster"] = 100*c["afdb_dark"].mean()
        s["lca_group_pct_of_clustered"] = {k: 100*v/c.height for k, v in c.group_by("afdb_lca_group").len().iter_rows()}
        s["lca_rank_pct_of_clustered"] = {str(k): 100*v/c.height for k, v in c.group_by("afdb_lca_rank").len().sort("len", descending=True).head(8).iter_rows()}
    s["pct_with_pfam_hit"] = 100*x["pfam_has_hit"].mean()
    s["pfam_n_regions_mean_among_hits"] = x.filter(pl.col("pfam_has_hit"))["pfam_n_regions"].mean()
    s["pct_pfam_hit_md5_matches_2026"] = 100*x.filter(pl.col("pfam_has_hit"))["pfam_seq_match_2026"].mean()
    s["pct_with_ted_or_pfam_modelled"] = 100*(m["ted_has_domain"] | m["pfam_has_hit"]).mean() if nm else None
    return s
summ = {"definitions": {
    "ted_has_domain": "protein has >=1 TED v1 domain (high or medium consensus) on its AFDB v4 model",
    "ted_class": "CATH class (alpha=1, beta=2, alphabeta=3) covering >=60% of TED-domain residues (unlabelled domains stay in denominator), else mixed; none = modelled but no TED domain; no_model = no AFDB v4 model",
    "afdb_in_cluster": "repId of file 5 is a cluster of AFDB v3 file 2 (non-singleton Foldseek cluster; includes AFDB50 non-representative members and clustered fragments)",
    "denominators": "pct_of_modelled_* use proteins with an AFDB v4 model; pct_of_afdb_entries_* use proteins present in AFDB v3 file 5"}}
summ["all"] = stats(a)
for p in ["euk", "prok"]: summ[p] = stats(a.filter(pl.col("panel") == p))
for sk in ["Eukaryota", "Bacteria", "Archaea"]: summ[sk] = stats(a.filter(pl.col("superkingdom") == sk))
summ["all_seq_match_only"] = stats(a.filter(pl.col("seq_match")))
summ["euk_seq_match_only"] = stats(a.filter((pl.col("panel") == "euk") & pl.col("seq_match")))
summ["prok_seq_match_only"] = stats(a.filter((pl.col("panel") == "prok") & pl.col("seq_match")))
json.dump(summ, open(R/"summary.json", "w"), indent=1)
print(json.dumps({k: summ[k] for k in ["euk", "prok"]}, indent=1)[:6000])
