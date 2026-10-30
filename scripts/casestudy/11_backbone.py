"""Print a backbone of the AA ML tree: collapse clades that are pure in one subfamily label; show UFBoot.
Rooted for display on the given outgroup label set (argv[2], comma-separated subfamilies) if it forms a clade,
else midpoint."""
import sys
import pandas as pd
from ete3 import Tree
R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
sf = sys.argv[1]
og = set(sys.argv[2].split(",")) if len(sys.argv) > 2 else set()
sel = pd.read_csv(f"{R}/results/casestudy/sample_members.tsv", sep="\t").set_index("seqid")
t = Tree(f"{R}/data/interim/casestudy/iqtree/{sf}_linsi_aa_ml.treefile", format=1)
lab = {l: sel.at[l, "subfamily"] + ("" if sel.at[l, "domain"] == "Eukaryota" else f"[{sel.at[l,'domain'][0]}]") for l in t.get_leaf_names()}
ogl = [l for l in t.get_leaf_names() if sel.at[l, "subfamily"] in og]
if ogl:
    anc = t.get_common_ancestor(ogl) if len(ogl) > 1 else t & ogl[0]
    if anc is t:
        t.set_outgroup(t & [l for l in t.get_leaf_names() if l not in ogl][0]); anc = t.get_common_ancestor(ogl)
    t.set_outgroup(anc)
else:
    t.set_outgroup(t.get_midpoint_outgroup())
t.ladderize()
for l in t.get_leaves():
    l.add_feature("lab", lab[l.name]); l.add_feature("k", 1)
for n in list(t.traverse("postorder")):
    if n.is_leaf():
        continue
    labs = {l.lab for l in n.get_leaves()}
    if len(labs) == 1:
        k = sum(l.k for l in n.get_leaves())
        bs = n.name
        for c in list(n.children):
            c.detach()
        n.add_feature("lab", labs.pop()); n.add_feature("k", k)
        n.name = f"{n.lab} x{k} (bs {bs})"
for l in t.get_leaves():
    if l.k == 1:
        l.name = l.lab
def show(n, d=0):
    if n.is_leaf():
        print("  " * d + f"- {n.name}  bl={n.dist:.2f}")
    else:
        print("  " * d + f"+ [{n.name}] bl={n.dist:.2f} n={len(n)}")
        for c in n.children: show(c, d + 1)
show(t)
