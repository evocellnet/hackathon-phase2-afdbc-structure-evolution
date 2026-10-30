"""X6 planning aid: pairwise p-distances (AA, 3Di) on the shared linsi columns, by pair class, to judge
saturation of pairwise ML distances for the WLS branch engine (deep prokaryote-eukaryote pairs)."""
import numpy as np, pandas as pd
R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
CS = f"{R}/data/interim/casestudy"
sel = pd.read_csv(f"{R}/results/casestudy/sample_members.tsv", sep="\t").set_index("seqid")
def read(p):
    d, k = {}, None
    for l in open(p):
        l = l.strip()
        if l.startswith(">"): k = l[1:]; d[k] = ""
        else: d[k] += l
    return d
rows = []
for sf in ["actin", "tubulin"]:
    A = read(f"{CS}/aln/{sf}_linsi_aa.fasta"); T = read(f"{CS}/aln/{sf}_linsi_3di.fasta")
    ks = list(A)
    a = np.array([list(A[k]) for k in ks]); t = np.array([list(T[k]) for k in ks])
    dom = np.array([sel.at[k, "domain"] for k in ks]); sub = np.array([sel.at[k, "subfamily"] for k in ks])
    euk = (dom == "Eukaryota") & (sub != "euk_FtsZ")
    for i in range(len(ks)):
        for j in range(i + 1, len(ks)):
            m = (a[i] != "-") & (a[j] != "-")
            n = m.sum()
            if n < 50: continue
            pa = (a[i][m] != a[j][m]).mean(); pt = (t[i][m] != t[j][m]).mean()
            cls = "euk-euk" if euk[i] and euk[j] else ("euk-nonEuk" if euk[i] or euk[j] else "nonEuk-nonEuk")
            rows.append((sf, cls, n, pa, pt))
df = pd.DataFrame(rows, columns=["superfamily", "pair_class", "n_cols", "p_aa", "p_3di"])
s = df.groupby(["superfamily", "pair_class"]).agg(n_pairs=("p_aa", "size"), median_p_aa=("p_aa", "median"),
    frac_p_aa_gt_0_75=("p_aa", lambda x: (x > 0.75).mean()), median_p_3di=("p_3di", "median"),
    frac_p_3di_gt_0_75=("p_3di", lambda x: (x > 0.75).mean())).round(3)
s.to_csv(f"{R}/results/casestudy/pairwise_pdist_by_class.tsv", sep="\t")
print(s.to_string())
