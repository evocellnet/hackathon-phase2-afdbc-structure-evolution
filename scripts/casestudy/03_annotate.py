"""X6 case study, step 1: superfamily membership + subfamily annotation.

Inputs
  data/interim/casestudy/hmm_chunks/c*.domtbl   hmmsearch --cut_ga of 26 Pfam HMMs (InterPro API, 2026-10-30)
                                                 vs all md5-matched panel proteins (euk 180 + prok 100)
  data/interim/casestudy/refs/refs.fasta        curated Swiss-Prot reference set (02_build_refs.py)
  data/processed/panel/protein_table*.parquet, species_panel*.tsv

Rules
  * actin superfamily: >= 1 GA hit to an actin-fold cytoskeletal family
    (Actin PF00022, MreB_Mbl PF06723, FtsA PF14450, StbA_N PF06406, ParM_N PF21523, ALP_N PF17989,
    MreB-like_C PF21522, Alp7A_like_C PF22128, Volactin PF25216, Salactin PF25229, Actin_micro PF17003)
    AND best cytoskeletal sequence score > best competitor score (HSP70 PF00012, PilM_2 PF11104).
  * tubulin superfamily: GA hit to the GTPase domain Tubulin PF00091 (Tubulin_2/Tubulin_3-only proteins,
    i.e. misato/DML1 and large bacterial multidomain proteins, are recorded but excluded)
    AND best tubulin score > Misat_Tub_SegII (misato) score.
  * quality: md5 match (already) and mean pLDDT >= 70; length within 0.5-2x the subfamily median.
  * subfamily: DIAMOND blastp (--ultra-sensitive) best hit (max bitscore) against the curated reference set,
    restricted to references of the same superfamily. Asgard (Promethearchaeota) and other archaeal/bacterial
    proteins whose best hit is a eukaryotic label are relabelled by taxonomy (e.g. 'Asgard_actin').
    Eukaryotic tubulins with best-hit identity < 40 % are labelled 'tubulin_divergent'
    (no zeta/eta tubulin is in Swiss-Prot).
"""
import glob, subprocess, os
import numpy as np, pandas as pd

R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
CS = f"{R}/data/interim/casestudy"
OUT = f"{R}/results/casestudy"
os.makedirs(OUT, exist_ok=True)

ACTIN = {"Actin", "MreB_Mbl", "FtsA", "StbA_N", "ParM_N", "ALP_N", "MreB-like_C", "Alp7A_like_C",
         "Volactin", "Salactin", "Actin_micro"}
ACTIN_COMP = {"HSP70", "PilM_2"}
TUB_CORE = {"Tubulin"}  # GTPase domain required; Tubulin_3 (misato) / Tubulin_2 (large bacterial proteins) alone do not qualify
TUB = TUB_CORE | {"Tubulin_C", "FtsZ_C", "CetZ_C", "TubZ_C", "TubZ_C_2", "TubZ_C_3", "TubZ-like_C", "FtsZ_C_2"}
TUB_COMP = {"Misat_Tub_SegII"}

# ---- parse hmmsearch ------------------------------------------------------
rows = []
for f in sorted(glob.glob(f"{CS}/hmm_chunks/c*.domtbl")):
    for line in open(f):
        if line.startswith("#"):
            continue
        x = line.split()
        rows.append((x[0], x[3], x[4], float(x[7]), int(x[2]), int(x[5]), float(x[13]), int(x[15]), int(x[16]),
                     int(x[17]), int(x[18])))
dom = pd.DataFrame(rows, columns=["target", "hmm", "pfam", "seq_score", "tlen", "qlen", "dom_score", "hmm_from",
                                  "hmm_to", "ali_from", "ali_to"])
dom.to_csv(f"{CS}/hmm_domtbl_all.tsv.gz", sep="\t", index=False)
seqsc = dom.groupby(["target", "hmm"]).seq_score.max().unstack(fill_value=0.0)
for h in ACTIN | ACTIN_COMP | TUB | TUB_COMP:
    if h not in seqsc:
        seqsc[h] = 0.0
best_act = seqsc[list(ACTIN)].max(axis=1)
best_actc = seqsc[list(ACTIN_COMP)].max(axis=1)
best_tub = seqsc[list(TUB)].max(axis=1)
best_tubcore = seqsc[list(TUB_CORE)].max(axis=1)
best_tubc = seqsc[list(TUB_COMP)].max(axis=1)
sf = pd.Series("", index=seqsc.index)
sf[(best_act > 0) & (best_act > best_actc)] = "actin"
tubmask = (best_tubcore > 0) & (best_tub > best_tubc)
both = tubmask & (sf == "actin")
print("proteins with both actin and tubulin hits:", int(both.sum()))
sf[tubmask & ~both] = "tubulin"
mem = pd.DataFrame({"target": seqsc.index, "superfamily": sf.values,
                    "best_cyto_score": np.where(sf.values == "tubulin", best_tub.values, best_act.values),
                    "competitor_score": np.where(sf.values == "tubulin", best_tubc.values, best_actc.values)})
# domain architecture string (ordered by envelope position) for members
arch = (dom.sort_values("ali_from").groupby("target").hmm.apply(lambda s: "+".join(dict.fromkeys(s))))
mem["pfam_arch"] = mem.target.map(arch)
excluded = mem[mem.superfamily == ""].copy()
mem = mem[mem.superfamily != ""].copy()
mem["acc"] = mem.target.str.split("|").str[0]
mem["leaf"] = mem.target.str.split("|").str[1]

# ---- panel metadata ------------------------------------------------------
pt = pd.concat([pd.read_parquet(f"{R}/data/processed/panel/{f}",
                                columns=["acc", "leaf", "db", "len_v4", "mean_plddt", "frac_plddt_ge70"])
                for f in ["protein_table.parquet", "protein_table_prok.parquet"]])
sp = pd.concat([pd.read_csv(f"{R}/data/processed/panel/{f}", sep="\t")
                for f in ["species_panel.tsv", "species_panel_prok.tsv"]])
sp = sp[["leaf", "taxid", "name", "superkingdom", "group", "lineage"]]
mem = mem.merge(pt, on=["acc", "leaf"], how="left").merge(sp, on="leaf", how="left")
mem["domain"] = mem.superkingdom.map({"Eukaryota": "Eukaryota", "Bacteria": "Bacteria", "Archaea": "Archaea"})
mem.loc[mem.group == "Promethearchaeota", "domain"] = "Asgard"
print(mem.groupby(["superfamily", "domain"]).size())

# ---- DIAMOND best hit to curated references ------------------------------
refs = f"{CS}/refs/refs.fasta"
lab = pd.read_csv(f"{CS}/refs/refs_labels.tsv", sep="\t")
seqs = {}
_st = set(mem.target)
with open(f"{CS}/members_all.fasta", "w") as fo:
    cur = None
    for line in open(f"{CS}/panel_all_v4match.fasta"):
        if line.startswith(">"):
            cur = line[1:].strip()
            if cur in _st:
                fo.write(line)
        elif cur in _st:
            fo.write(line)
            seqs[cur] = seqs.get(cur, "") + line.strip()
subprocess.run(f"diamond makedb --in {refs} -d {CS}/refs/refs --quiet", shell=True, check=True)
subprocess.run(f"diamond blastp -q {CS}/members_all.fasta -d {CS}/refs/refs -o {CS}/members_vs_refs.tsv "
               f"--ultra-sensitive -k 25 -e 1e-5 --threads 8 --quiet "
               f"--outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore qcovhsp", shell=True, check=True)
bh = pd.read_csv(f"{CS}/members_vs_refs.tsv", sep="\t", header=None,
                 names=["target", "ref", "pident", "alen", "qlen", "slen", "evalue", "bits", "qcov"])
bh["ref_label"] = bh.ref.str.split("|").str[1]
bh = bh.merge(lab[["label", "superfamily"]].drop_duplicates("label").rename(
    columns={"label": "ref_label", "superfamily": "ref_sf"}), on="ref_label")
bh = bh.merge(mem[["target", "superfamily"]], on="target")
bh = bh[bh.ref_sf == bh.superfamily].sort_values(["target", "bits"], ascending=[True, False])
top = bh.groupby("target").head(1).set_index("target")
# second-best distinct label, to report ambiguity
sec = (bh[bh.ref_label != bh.target.map(top.ref_label)].groupby("target").head(1).set_index("target"))
mem["best_ref"] = mem.target.map(top.ref.str.split("|").str[0])
mem["best_label"] = mem.target.map(top.ref_label)
mem["best_pident"] = mem.target.map(top.pident)
mem["best_bits"] = mem.target.map(top.bits)
mem["second_label"] = mem.target.map(sec.ref_label)
mem["second_bits"] = mem.target.map(sec.bits)

EUK_LABELS = lambda l: isinstance(l, str) and (l.startswith("tubulin") or l.startswith("Arp") or l.startswith("actin"))


PFAM_FALLBACK = [("Salactin", "salactin"), ("Volactin", "volactin"), ("ParM_N", "ALP_bacterial"),
                 ("StbA_N", "ALP_bacterial"), ("ALP_N", "ALP_bacterial"), ("FtsA", "FtsA"), ("MreB_Mbl", "MreB"),
                 ("Actin", "actin_divergent")]


def subfam(r):
    l = r.best_label
    if not isinstance(l, str):
        if r.superfamily == "actin":
            for h, lab_ in PFAM_FALLBACK:
                if h in str(r.pfam_arch).split("+"):
                    l = lab_
                    break
        if not isinstance(l, str):
            return f"{r.superfamily}_nohit"
        if r.domain == "Eukaryota" and l in ("MreB", "FtsA"):
            return f"euk_{l}"
        if r.domain == "Asgard":
            return "Asgard_actin"
        return l
    if r.domain == "Bacteria" and l == "Ta0583_archaeal_actin":
        return "ALP_bacterial"
    if r.domain == "Archaea" and r.group == "Nitrososphaerota" and l.startswith("tubulin"):
        return "artubulin"
    if r.domain == "Asgard":
        if r.superfamily == "actin":
            return "Asgard_actin" if EUK_LABELS(l) or l == "crenactin" else f"Asgard_{l}"
        return "Asgard_tubulin" if EUK_LABELS(l) else f"Asgard_{l}"
    if r.domain != "Eukaryota" and EUK_LABELS(l):
        return f"prok_{r.superfamily}-like({l})"
    if r.domain == "Eukaryota" and l.startswith("tubulin") and r.best_pident < 40:
        return "tubulin_divergent"
    if r.domain == "Eukaryota" and l in ("BtubA", "BtubB"):
        return "tubulin_divergent"
    if r.domain == "Eukaryota" and not EUK_LABELS(l):
        return f"euk_{l}"  # e.g. organellar FtsZ in plants/algae
    return l


mem["subfamily"] = mem.apply(subfam, axis=1)
mem["plddt_ok"] = mem.mean_plddt >= 70
med = mem[mem.plddt_ok].groupby("subfamily").len_v4.median()
mem["len_ok"] = (mem.len_v4 >= 0.5 * mem.subfamily.map(med)) & (mem.len_v4 <= 2 * mem.subfamily.map(med))
mem["pass"] = mem.plddt_ok & mem.len_ok
cols = ["acc", "leaf", "taxid", "name", "domain", "group", "superfamily", "subfamily", "best_label", "best_ref",
        "best_pident", "best_bits", "second_label", "second_bits", "pfam_arch", "best_cyto_score",
        "competitor_score", "db", "len_v4", "mean_plddt", "frac_plddt_ge70", "plddt_ok", "len_ok", "pass", "target"]
mem = mem[cols].sort_values(["superfamily", "domain", "subfamily", "group", "name"])
mem.to_csv(f"{OUT}/membership_all.tsv", sep="\t", index=False)
excluded.to_csv(f"{CS}/hmm_hits_not_members.tsv", sep="\t", index=False)
summ = (mem.groupby(["superfamily", "domain", "subfamily"])
        .agg(n_hits=("acc", "size"), n_pass=("pass", "sum"), n_species=("leaf", "nunique"))
        .reset_index())
summ.to_csv(f"{OUT}/membership_summary.tsv", sep="\t", index=False)
pd.set_option("display.width", 200); pd.set_option("display.max_rows", 300)
print(summ.to_string())
