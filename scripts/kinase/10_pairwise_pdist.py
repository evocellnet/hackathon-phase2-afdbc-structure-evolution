"""X7 saturation check: pairwise AA and 3Di p-distances (shared non-gap columns of the filtered alignment) summarised
by pair class (tip categories / ePK groups). Usage: python 10_pairwise_pdist.py <src>"""
import sys
import numpy as np, pandas as pd
R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
K = f"{R}/data/interim/kinase"
src = sys.argv[1] if len(sys.argv) > 1 else "linsi"
sel = pd.read_csv(f"{R}/results/kinase/sample_members.tsv", sep="\t").set_index("seqid")


def rf(p):
    d, k = {}, None
    for line in open(p):
        line = line.strip()
        if line.startswith(">"):
            k = line[1:]; d[k] = []
        elif k:
            d[k].append(line)
    return {k: "".join(v) for k, v in d.items()}


A, T = rf(f"{K}/aln/{src}_aa.fasta"), rf(f"{K}/aln/{src}_3di.fasta")
keys = list(A)
Ma = np.array([list(A[k]) for k in keys]); Mt = np.array([list(T[k]) for k in keys])
ok = Ma != "-"


def cat(c):
    if c.startswith("ePK_"):
        return "euk_ePK:" + c[4:]
    if c in ("prok_ePK", "Asgard_ePK"):
        return "prok_ePK"
    return ("euk_ELK:" if c.endswith("_euk") else "prok_ELK:") + c.rsplit("_", 1)[0]


C = [cat(sel.at[k, "klass"]) for k in keys]
rows = []
for i in range(len(keys)):
    both = ok[i] & ok[i + 1:]
    n = both.sum(1)
    pa = ((Ma[i] != Ma[i + 1:]) & both).sum(1) / np.maximum(n, 1)
    pt = ((Mt[i] != Mt[i + 1:]) & both).sum(1) / np.maximum(n, 1)
    for j in range(len(n)):
        a, b = sorted([C[i], C[i + 1 + j]])
        ca, cb = a.split(":")[0], b.split(":")[0]
        if ca == cb == "euk_ePK":
            pc = "euk_ePK same group" if a == b else "euk_ePK different group"
        elif ca == cb and ca != "prok_ePK":
            pc = f"{ca} same family" if a == b else f"{ca} different family"
        else:
            pc = f"{ca} x {cb}"
        rows.append((pc, n[j], pa[j], pt[j]))
d = pd.DataFrame(rows, columns=["pair_class", "n_cols", "p_aa", "p_3di"])
d = d[d.n_cols >= 50]
s = d.groupby("pair_class").agg(n_pairs=("p_aa", "size"), median_p_aa=("p_aa", "median"),
                                frac_p_aa_gt_0p75=("p_aa", lambda x: (x > 0.75).mean()),
                                median_p_3di=("p_3di", "median"),
                                frac_p_3di_gt_0p75=("p_3di", lambda x: (x > 0.75).mean())).round(3)
s.to_csv(f"{R}/results/kinase/pairwise_pdist_by_class_{src}.tsv", sep="\t")
print(s.to_string())
