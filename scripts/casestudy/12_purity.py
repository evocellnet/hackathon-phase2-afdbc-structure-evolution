"""Per-subfamily coherence on the AA ML tree: largest clade (bipartition side) made only of that label, as a fraction
of the label's sampled members, and the UFBoot of that clade. Complements clade_checks_*.tsv."""
import pandas as pd
from ete3 import Tree
R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
sel = pd.read_csv(f"{R}/results/casestudy/sample_members.tsv", sep="\t").set_index("seqid")
rows = []
for sf in ["actin", "tubulin"]:
    t = Tree(f"{R}/data/interim/casestudy/iqtree/{sf}_linsi_aa_ml.treefile", format=1)
    leaves = set(t.get_leaf_names())
    sides = []
    for n in t.traverse():
        if n.is_root(): continue
        s = set(n.get_leaf_names())
        sup = n.name if not n.is_leaf() else ""
        sides += [(frozenset(s), sup), (frozenset(leaves - s), sup)]
    sub = sel.loc[list(leaves), "subfamily"]
    for lab, mem in sub.groupby(sub):
        m = set(mem.index)
        best, bsup = 0, ""
        for s, sup in sides:
            if s <= m and len(s) > best:
                best, bsup = len(s), sup
        dom = sel.loc[list(m), "domain"].iloc[0]
        rows.append(dict(superfamily=sf, subfamily=lab, domain=dom, n=len(m), largest_pure_clade=best,
                         frac=round(best / len(m), 2), ufboot=bsup))
df = pd.DataFrame(rows)
df.to_csv(f"{R}/results/casestudy/subfamily_coherence.tsv", sep="\t", index=False)
print(df.to_string())
