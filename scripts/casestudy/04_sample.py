"""X6 step 2: downsample each superfamily to <= 400 sequences (seed 20261030).

* Only proteins with pass == True (md5 match, mean pLDDT >= 70, length 0.5-2x subfamily median).
* Keep ALL prokaryotic (Bacteria, Archaea) and Asgard members.
* Eukaryotic members of likely-contaminant labels (euk_MreB, euk_FtsA: bacterial-type MreB/FtsA in animal
  genomes, e.g. Glossina endosymbionts) are not sampled.
* Eukaryotes fill the remaining budget by stratified round-robin over cells (subfamily x supergroup):
  each round draws one random member per non-empty cell, preferring species not yet drawn for that subfamily;
  the last round is drawn in random cell order until the budget is reached. Lineage-restricted orphan
  labels are capped (CAP).
"""
import numpy as np, pandas as pd

R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
CS = f"{R}/data/interim/casestudy"
OUT = f"{R}/results/casestudy"
SEED = 20261030
BUDGET = 400
EXCL = {"euk_MreB", "euk_FtsA"}
CAP = {"actin_POTE": 3, "actin-like_metazoan_orphan": 6, "actin-like_other": 6, "actin_divergent": 2,
       "euk_FtsZ": 24, "tubulin_divergent": 16}

m = pd.read_csv(f"{OUT}/membership_all.tsv", sep="\t")
m = m[m["pass"]].copy()
rng = np.random.default_rng(SEED)
sel_all = []
for sf in ["actin", "tubulin"]:
    d = m[m.superfamily == sf]
    keep = d[d.domain != "Eukaryota"]
    e = d[(d.domain == "Eukaryota") & ~d.subfamily.isin(EXCL)].copy()
    e = e.sample(frac=1.0, random_state=SEED)  # deterministic shuffle
    budget = BUDGET - len(keep)
    cells = {k: list(g.index) for k, g in e.groupby(["subfamily", "group"], sort=True)}
    chosen, used_species, n_sub = [], {}, {}
    while len(chosen) < budget and any(cells.values()):
        order = list(cells.keys())
        rng.shuffle(order)
        for c in order:
            if len(chosen) >= budget:
                break
            sub = c[0]
            if not cells[c] or n_sub.get(sub, 0) >= CAP.get(sub, 10 ** 9):
                continue
            idx = cells[c]
            us = used_species.setdefault(sub, set())
            pick = next((i for i in idx if e.at[i, "leaf"] not in us), idx[0])
            idx.remove(pick)
            chosen.append(pick)
            us.add(e.at[pick, "leaf"])
            n_sub[sub] = n_sub.get(sub, 0) + 1
        if all((not v) or n_sub.get(k[0], 0) >= CAP.get(k[0], 10 ** 9) for k, v in cells.items()):
            break
    s = pd.concat([keep, e.loc[chosen]])
    s["sampled_as"] = np.where(s.domain == "Eukaryota", "stratified", "all_kept")
    sel_all.append(s)
    print(sf, len(s), "(non-euk", len(keep), ")")
    print(s.groupby(["domain", "subfamily"]).size().to_string())
sel = pd.concat(sel_all)
sel["seqid"] = sel.acc + "_" + sel.taxid.astype(str)
sel.to_csv(f"{OUT}/sample_members.tsv", sep="\t", index=False)
seqs, cur = {}, None
want = set(sel.target)
for line in open(f"{CS}/members_all.fasta"):
    if line.startswith(">"):
        cur = line[1:].strip()
    elif cur in want:
        seqs[cur] = seqs.get(cur, "") + line.strip()
for sf in ["actin", "tubulin"]:
    with open(f"{CS}/{sf}_sample.fasta", "w") as fo:
        for _, r in sel[sel.superfamily == sf].iterrows():
            fo.write(f">{r.seqid}\n{seqs[r.target]}\n")
print(sel.groupby(["superfamily", "group"]).size().unstack(0).fillna(0).astype(int).to_string())
