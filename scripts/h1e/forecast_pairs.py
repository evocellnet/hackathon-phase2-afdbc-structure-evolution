#!/usr/bin/env python
"""Forecast of ortholog pairs with qualifying PDB chains on both sides (symbol-matched vertebrate pairs; proteome fractions for all taxa).
Writes data/processed/h1e/pair_forecast.tsv and proteome_fractions.tsv"""
import polars as pl, itertools
D = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution/data"
TAXG = {9606: "human", 10090: "mouse", 10116: "rat", 9913: "bovine", 9031: "chicken", 7955: "zebrafish"}
rp = pl.read_csv(f"{D}/raw/proteomes/meta/refprot_meta.tsv", separator="\t").select(pl.col("Organism Id").alias("taxid"), pl.col("Protein count").alias("proteome_size"), "Organism")
# denominators: entries per taxon in Foldseek afdb_proteome (Swiss-Prot+TrEMBL, so human ~42k not 20k); bovine/chicken are absent there -> refprot_meta protein count
am = pl.read_csv(f"{D}/raw/foldseek_db/afdb_proteome/afdb_proteome_mapping", separator="\t", has_header=False, new_columns=["key", "taxid"]).group_by("taxid").len()
am = {r[0]: r[1] for r in am.iter_rows()}
T = {9606: "human", 10090: "mouse", 10116: "rat", 559292: "yeast", 3702: "athaliana", 7227: "dmel", 6239: "celegans", 7955: "zebrafish", 83333: "ecoli_k12", 9031: "chicken", 9913: "bovine"}
out = {}
for key in ("pdb_chains.tsv", "pdb_chains_post_af2.tsv"):
    c = pl.read_csv(f"{D}/processed/h1e/{key}", separator="\t", infer_schema_length=100000, schema_overrides={"taxid": pl.Int64}).filter(pl.col("qualifies") == True)
    out[key] = c
frac = []
for t, n in T.items():
    ps = am.get(t)
    if ps is None:
        ps = rp.filter(pl.col("taxid") == t)["proteome_size"]; ps = int(ps[0]) if len(ps) else None
    r = [t, n, ps]
    for key in out:
        k = out[key].filter(pl.col("taxid") == t).height; r += [k, round(k / ps, 4) if ps else None]
    frac.append(r)
pl.DataFrame(frac, schema=["taxid", "name", "proteome_size", "n_chains", "frac", "n_postAF2", "frac_postAF2"], orient="row").write_csv(f"{D}/processed/h1e/proteome_fractions.tsv", separator="\t")
G = {}
for t in TAXG:
    g = pl.read_csv(f"{D}/raw/pdb_meta/uniprot_3d_genes_{t}.tsv", separator="\t").select(pl.col("Entry").alias("uniprot"), pl.col("Gene Names (primary)").alias("gene").str.to_uppercase())
    G[t] = g.filter(pl.col("gene").is_not_null())
rows = []
for a, b in itertools.combinations(TAXG, 2):
    r = [TAXG[a], TAXG[b]]
    for key, c in out.items():
        sets = []
        for t in (a, b):
            sets.append(set(c.filter(pl.col("taxid") == t).join(G[t], on="uniprot")["gene"].to_list()))
        r += [len(sets[0]), len(sets[1]), len(sets[0] & sets[1])]
    rows.append(r)
pl.DataFrame(rows, schema=["A", "B", "nA", "nB", "pairs_both", "nA_post", "nB_post", "pairs_both_postAF2"], orient="row").write_csv(f"{D}/processed/h1e/pair_forecast.tsv", separator="\t")
pl.Config.set_tbl_rows(40); pl.Config.set_tbl_cols(12)
print(pl.read_csv(f"{D}/processed/h1e/proteome_fractions.tsv", separator="\t")); print(pl.read_csv(f"{D}/processed/h1e/pair_forecast.tsv", separator="\t"))
