#!/usr/bin/env python
"""Build data/processed/h1e/pdb_chains.tsv (+ pdb_chains_post_af2.tsv): best experimental chain per UniProt accession.

Inputs (data/raw/pdb_meta): pdb_chain_uniprot.tsv.gz (SIFTS), pdb_entry_meta.tsv.gz (RCSB GraphQL),
uniprot_3d_<taxid>.tsv (UniProt length + taxon membership for the 11 taxa), Foldseek afdb_proteome index (lengths for 48 panel proteomes).
Accessions outside these sets keep coverage NA (ranked by n_mapped).
Filters: X-ray, or EM with res <= 3.0; resolution <= 3.0 A; protein UniProt canonical accessions (isoforms dropped).
Ranking: UniProt coverage (desc; = n_mapped residues, same accession) then resolution (asc) then release date (desc).
Coverage = min(merged SIFTS segment residues, RCSB modeled residue count) of the chain / UniProt length (observed-residue correction;
the uncorrected SIFTS construct coverage is kept as coverage_sifts_construct). 'qualifies' = coverage >= 0.70.
Two selections: (a) best overall chain (post_af2_cutoff flag = release > 2018-04-30);
(b) best chain among entries released after 2018-04-30 only (DESIGN 11.10)."""
import polars as pl, os
R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution/data/raw/pdb_meta"
O = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution/data/processed/h1e"
CUT = "2018-04-30"; MINCOV = 0.70; MAXRES = 3.0
TAX = {9606: "human", 10090: "mouse", 10116: "rat", 559292: "yeast", 3702: "athaliana", 7227: "dmel", 6239: "celegans",
       7955: "zebrafish", 83333: "ecoli_k12", 9031: "chicken", 9913: "bovine"}
s = pl.read_csv(f"{R}/pdb_chain_uniprot.tsv.gz", separator="\t", skip_rows=1, infer_schema_length=0)
s = s.rename({"PDB": "pdb_id", "CHAIN": "chain", "SP_PRIMARY": "uniprot"}).with_columns(
    pl.col("SP_BEG").cast(pl.Int64), pl.col("SP_END").cast(pl.Int64))
s = s.filter(~pl.col("uniprot").str.contains("-"))
# merged covered residues per (pdb, chain, uniprot): vectorised interval union (sort by start; new residues = end - max(start-1, running max end))
k = ["pdb_id", "chain", "uniprot"]
s = s.sort(k + ["SP_BEG"]).with_columns(pl.col("SP_END").cum_max().over(k).shift(1).over(k).fill_null(0).alias("prev_end"))
s = s.with_columns(pl.max_horizontal(pl.lit(0), pl.col("SP_END") - pl.max_horizontal(pl.col("SP_BEG") - 1, pl.col("prev_end"))).alias("new_res"))
c = s.group_by(k).agg(pl.col("SP_BEG").min().alias("sp_beg"), pl.col("SP_END").max().alias("sp_end"), pl.col("new_res").sum().alias("n_mapped"),
                      pl.len().alias("n_segments"))
m = pl.read_csv(f"{R}/pdb_entry_meta.tsv.gz", separator="\t", infer_schema_length=0).with_columns(
    pl.col("resolution").cast(pl.Float64, strict=False)).rename({"rcsb_method": "method", "initial_release_date": "release_date"})
c = c.join(m.select("pdb_id", "method", "resolution", "release_date"), on="pdb_id", how="inner")
c = c.filter(pl.col("method").is_in(["X-ray", "EM"]) & (pl.col("resolution") <= MAXRES))
# SIFTS segments span SEQRES incl. residues without coordinates -> cap by the chain's modeled residue count (RCSB; upper bound of
# observed mapped residues since it also counts tags/linkers). n_sifts = SIFTS construct residues, n_mapped = min(n_sifts, modeled).
md = pl.read_csv(f"{R}/pdb_chain_modeled.tsv.gz", separator="\t", schema_overrides={"modeled_residue_count": pl.Int64, "entity_seq_length": pl.Int64, "chain": pl.String}).unique(["pdb_id", "chain"])
c = c.join(md.select("pdb_id", "chain", "modeled_residue_count"), on=["pdb_id", "chain"], how="left")
c = c.rename({"n_mapped": "n_sifts"}).with_columns(pl.min_horizontal("n_sifts", "modeled_residue_count").fill_null(pl.col("n_sifts")).alias("n_mapped"))
# UniProt lengths: (i) per-taxon UniProt REST tables for the 11 H1e taxa (structure_3d:true); (ii) Foldseek afdb_proteome (48 panel
# proteomes; seq length = .index col3 - 2, accession from .lookup, taxid from _mapping) for the remaining panel species.
ul = [pl.read_csv(f"{R}/uniprot_3d_{t}.tsv", separator="\t", schema_overrides={"Length": pl.Int64, "Organism (ID)": pl.Int64})
      .rename({"Entry": "uniprot", "Length": "uniprot_length", "Organism (ID)": "organism_id", "Reviewed": "reviewed"}) for t in TAX]
FS = f"{R}/../foldseek_db/afdb_proteome/afdb_proteome"
fl = pl.read_csv(f"{FS}.lookup", separator="\t", has_header=False, new_columns=["key", "afid", "f"]).with_columns(
    pl.col("afid").str.extract(r"^AF-(.+)-F\d+-model", 1).alias("uniprot"))
fi = pl.read_csv(f"{FS}.index", separator="\t", has_header=False, new_columns=["key", "off", "len"]).with_columns((pl.col("len") - 2).alias("uniprot_length"))
fm = pl.read_csv(f"{FS}_mapping", separator="\t", has_header=False, new_columns=["key", "organism_id"])
fs_ = fl.join(fi.select("key", "uniprot_length"), on="key").join(fm, on="key").select("uniprot", "uniprot_length", "organism_id").with_columns(pl.lit(None, dtype=pl.String).alias("reviewed"))
ul = pl.concat([x.select("uniprot", "uniprot_length", "organism_id", "reviewed") for x in ul] + [fs_]).unique("uniprot", keep="first")
c = c.join(ul, on="uniprot", how="left")
mem = []
for t in TAX:
    f = f"{R}/uniprot_3d_{t}.tsv"
    if os.path.exists(f): mem.append(pl.read_csv(f, separator="\t").select(pl.col("Entry").alias("uniprot")).with_columns(pl.lit(t).alias("taxid")))
mem = pl.concat(mem).unique("uniprot")
c = c.join(mem, on="uniprot", how="left")
c = c.with_columns((pl.col("n_mapped") / pl.col("uniprot_length")).alias("coverage"), (pl.col("n_sifts") / pl.col("uniprot_length")).alias("coverage_sifts_construct"),
                   (pl.col("release_date") > CUT).alias("post_af2_cutoff"))
c = c.with_columns((pl.col("coverage") >= MINCOV).alias("qualifies"))
print("candidate chain rows", c.height, "accessions", c["uniprot"].n_unique(), "no UniProt length", c.filter(pl.col("uniprot_length").is_null())["uniprot"].n_unique())
def best(df):
    return (df.filter(pl.col("coverage").is_null() | pl.col("qualifies"))
              .sort(["uniprot", "n_mapped", "resolution", "release_date"], descending=[False, True, False, True])
              .unique("uniprot", keep="first", maintain_order=True))
cols = ["uniprot", "taxid", "organism_id", "reviewed", "pdb_id", "chain", "method", "resolution", "release_date", "post_af2_cutoff",
        "sp_beg", "sp_end", "n_sifts", "modeled_residue_count", "n_mapped", "n_segments", "uniprot_length", "coverage", "coverage_sifts_construct", "qualifies"]
a = best(c).select(cols).sort("uniprot"); b = best(c.filter(pl.col("post_af2_cutoff"))).select(cols).sort("uniprot")
os.makedirs(O, exist_ok=True)
a.write_csv(f"{O}/pdb_chains.tsv", separator="\t", float_precision=4); b.write_csv(f"{O}/pdb_chains_post_af2.tsv", separator="\t", float_precision=4)
print("pdb_chains.tsv", a.height, "post_af2 table", b.height)
print(a.filter(pl.col("qualifies").is_null()).height, "rows lack UniProt length (coverage NA)")
# taxon counts
out = []
for t, nm in TAX.items():
    r = [t, nm]
    for d in (a, b):
        q = d.filter((pl.col("taxid") == t) & pl.col("qualifies"))
        r += [q.height, q.filter(pl.col("method") == "X-ray").height]
    out.append(r)
T = pl.DataFrame(out, schema=["taxid", "name", "n_any", "n_xray", "n_postAF2_any", "n_postAF2_xray"], orient="row")
T.write_csv(f"{O}/taxon_counts.tsv", separator="\t"); print(T)
