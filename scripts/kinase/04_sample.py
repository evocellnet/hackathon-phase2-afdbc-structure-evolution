"""X7 step 4: downsample the PKL tree set to <= 500 kinase domains (seed 20261030), stratified by group.

Only 'pass' domains (HMM coverage >= 0.7, domain pLDDT >= 70, domain length <= 2x family median).
Strata (klass from 03_annotate.py) and caps:
  * prokaryotic ePK-like (Pkinase/PK_Tyr best hit in Bacteria/Archaea): <= 80, round-robin over (domain, phylum)
  * Asgard ePK-like: all
  * ELK/atypical classes RIO1/ABC1/APH/Choline_kinase/Kdo x {euk, bac, arc}: <= 12 each, round-robin over group
  * eukaryotic best-hit-TK outside Holozoa (ePK_TKlike_nonholozoan): <= 10
  * eukaryotic ePKs: the remaining budget split equally over AGC, CAMK, CK1, CMGC, STE, TK, TKL, Other
    (ePK_unassigned not sampled); within a group, round-robin over supergroups, one random member per cell per round,
    preferring species not yet drawn for that group.
Outputs results/kinase/sample_members.tsv, data/interim/kinase/pkl_sample.fasta (domain sequences; id acc_taxid)
"""
import numpy as np, pandas as pd

R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
K = f"{R}/data/interim/kinase"
OUT = f"{R}/results/kinase"
SEED, BUDGET = 20261030, 500
EPK_GROUPS = ["AGC", "CAMK", "CK1", "CMGC", "STE", "TK", "TKL", "Other"]
ELK = ["RIO1", "ABC1", "APH", "Choline_kinase", "Kdo"]
rng = np.random.default_rng(SEED)
m = pd.read_csv(f"{OUT}/membership_all.tsv.gz", sep="\t")
m = m[m["pass"]].sample(frac=1.0, random_state=SEED)


def rr(d, cap, cellcol="group"):
    """round-robin over cells (one random member per cell per round, preferring unused species)"""
    cells = {k: list(g.index) for k, g in d.groupby(cellcol, sort=True)}
    chosen, used = [], set()
    while len(chosen) < cap and any(cells.values()):
        order = sorted(cells)
        rng.shuffle(order)
        for c in order:
            if len(chosen) >= cap or not cells[c]:
                continue
            idx = cells[c]
            pick = next((i for i in idx if d.at[i, "leaf"] not in used), idx[0])
            idx.remove(pick); chosen.append(pick); used.add(d.at[pick, "leaf"])
    return chosen


parts = []
p = m[m.klass == "prok_ePK"].copy(); p["cell"] = p.domain + ":" + p.group
parts.append(m.loc[rr(p, 80, "cell")])
parts.append(m[m.klass == "Asgard_ePK"])
for f in ELK:
    for d in ("euk", "bac", "arc", "asg"):
        x = m[m.klass == f"{f}_{d}"]
        if len(x):
            parts.append(m.loc[rr(x, 12)])
parts.append(m.loc[rr(m[m.klass == "ePK_TKlike_nonholozoan"], 10)])
n_fixed = sum(len(x) for x in parts)
per = (BUDGET - n_fixed) // len(EPK_GROUPS)
rest = (BUDGET - n_fixed) - per * len(EPK_GROUPS)
for i, g in enumerate(EPK_GROUPS):
    x = m[m.klass == f"ePK_{g}"]
    parts.append(m.loc[rr(x, per + (1 if i < rest else 0))])
s = pd.concat(parts)
s["seqid"] = s.acc + "_" + s.taxid.astype(str)
assert s.seqid.is_unique and len(s) <= BUDGET
s.to_csv(f"{OUT}/sample_members.tsv", sep="\t", index=False)
seqs, cur = {}, None
for line in open(f"{K}/pkl_domains.fasta"):
    if line.startswith(">"):
        cur = line[1:].strip()
    else:
        seqs[cur] = seqs.get(cur, "") + line.strip()
with open(f"{K}/pkl_sample.fasta", "w") as fo:
    for r in s.itertuples():
        fo.write(f">{r.seqid}\n{seqs[r.target]}\n")
pd.set_option("display.width", 200)
print(len(s)); print(s.groupby(["klass"]).agg(n=("acc", "size"), species=("leaf", "nunique"),
                                              groups=("group", "nunique")).to_string())
print(s[s.domain == "Eukaryota"].groupby(["group"]).size().to_string())
