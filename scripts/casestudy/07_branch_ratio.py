"""X6 step 4 (first look): per-branch b_3Di / b_AA on the fixed AA ML topology.

Branch lengths: IQ-TREE3 -te <ML tree> with LG+G4 on the AA columns and Q.3Di.AF+G4 on the 3Di letters
projected onto the SAME columns (06_alignments.py). Branches are matched by bipartition.
rho_e = (b_3Di,e / b_AA,e) / (sum b_3Di / sum b_AA)  -- the 3Di and AA scales differ, so only the
relative ratio is interpretable. Internal branches are 'supported' if UFBoot >= 95 (DESIGN 11.8).
Focal branches: the stem of the eukaryotic clade (actin+Arps; eukaryotic tubulins) and its relation to
Asgard / prokaryotic homologues. Descriptive only, no tests (11.8).
Usage: python 07_branch_ratio.py <sf> <src> <tag>
"""
import sys, json
import numpy as np, pandas as pd
from ete3 import Tree

R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
CS = f"{R}/data/interim/casestudy"
OUT = f"{R}/results/casestudy"
sf, src = sys.argv[1], sys.argv[2]
tag = sys.argv[3] if len(sys.argv) > 3 else "binary"
T = f"{CS}/iqtree"
ml = Tree(f"{T}/{sf}_{src}_aa_ml.treefile", format=1)
aa = Tree(f"{T}/{sf}_{src}_aa_fixed_{tag}.treefile", format=1)
td = Tree(f"{T}/{sf}_{src}_3di_fixed_{tag}.treefile", format=1)
sel = pd.read_csv(f"{OUT}/sample_members.tsv", sep="\t").set_index("seqid")
leaves = sorted(ml.get_leaf_names())
ref = leaves[0]
ALL = frozenset(leaves)


def splits(t, with_support=False):
    out = {}
    for n in t.traverse():
        if n.is_root():
            continue
        s = frozenset(n.get_leaf_names())
        if ref in s:
            s = ALL - s
        sup = np.nan
        if with_support and not n.is_leaf():
            try:
                sup = float(str(n.name).split("/")[-1]) if n.name else np.nan
            except ValueError:
                sup = np.nan
        out[s] = (n.dist, sup, n.is_leaf())
    return out


if src == "linsi":  # write the §11.8 reference topology: collapse internal branches with UFBoot < 95
    col = ml.copy()
    for n in list(col.traverse("postorder")):
        if n.is_leaf() or n.is_root():
            continue
        try:
            sup = float(str(n.name).split("/")[-1])
        except ValueError:
            sup = 0.0
        if sup < 95:
            n.delete(preserve_branch_length=False)
    col.write(outfile=f"{OUT}/{sf}_reference_topology_ufb95.nwk", format=1)
    ml.write(outfile=f"{OUT}/{sf}_aa_ml_ufboot.nwk", format=1)

S_ml, S_aa, S_td = splits(ml, True), splits(aa), splits(td)
# merge root-adjacent duplicate splits (unrooted: the two root edges are one bipartition)
rows = []
for s in set(S_aa) | set(S_td):
    if len(s) == 0 or len(s) == len(ALL):
        continue
    b_aa = S_aa.get(s, (np.nan,))[0]
    b_td = S_td.get(s, (np.nan,))[0]
    sup = S_ml.get(s, (np.nan, np.nan, False))[1]
    leaf = len(s) == 1 or len(s) == len(ALL) - 1
    side = s if len(s) <= len(ALL) / 2 else ALL - s
    dom = sel.loc[list(side), "domain"].value_counts().to_dict()
    sub = sel.loc[list(side), "subfamily"].value_counts().to_dict()
    rows.append(dict(split_size=len(side), terminal=leaf, ufboot=sup, b_aa=b_aa, b_3di=b_td,
                     tip=list(side)[0] if len(side) == 1 else "",
                     domains=json.dumps(dom), subfamilies=json.dumps(dict(sorted(sub.items(), key=lambda x: -x[1])[:6])),
                     _side=side))
df = pd.DataFrame(rows)
in_both = df.b_aa.notna() & df.b_3di.notna()
assert in_both.all(), "topologies differ between fixed-tree runs"
scale = df.b_3di.sum() / df.b_aa.sum()
df["ratio"] = df.b_3di / df.b_aa
df["rho"] = df.ratio / scale
df["supported"] = df.terminal | (df.ufboot >= 95)

# ---- classify branches ------------------------------------------------------------------------
euk_excl = {"euk_FtsZ"}
dom_of = sel.domain.to_dict(); sub_of = sel.subfamily.to_dict()
EUK = frozenset(k for k in leaves if dom_of[k] == "Eukaryota" and sub_of[k] not in euk_excl)
ASG = frozenset(k for k in leaves if dom_of[k] == "Asgard")


def cls(side):
    d = {dom_of[k] for k in side}
    other = ALL - side
    if side == EUK or other == EUK:
        return "EUK_STEM"
    if ASG and (side == EUK | ASG or other == EUK | ASG):
        return "EUK+ASGARD_STEM"
    sd = side if side <= EUK or side.isdisjoint(EUK) else other
    if side <= EUK or other <= EUK:
        return "within_eukaryote"
    if side.isdisjoint(EUK) or other.isdisjoint(EUK):
        return "within_prok_or_mixed"
    return "mixed_split"


df["class"] = [cls(s) for s in df._side]
# subfamily clade stems (exact monophyly of a eukaryotic subfamily)
lab = {}
for sub in sorted(set(sub_of.values())):
    members = frozenset(k for k in leaves if sub_of[k] == sub)
    if len(members) < 2:
        continue
    for i, s in enumerate(df._side):
        if s == members or (ALL - s) == members:
            lab[i] = f"stem:{sub}"
df["clade_label"] = [lab.get(i, "") for i in range(len(df))]


def min_clade(members):
    """smallest side of any bipartition containing all members (unrooted 'clade' check)."""
    best = None
    for s in df._side:
        for side in (s, ALL - s):
            if members <= side and (best is None or len(side) < len(best)):
                best = side
    return best


checks = []
groups = {f"sub:{s}": frozenset(k for k in leaves if sub_of[k] == s) for s in set(sub_of.values())}
if sf == "actin":
    groups["MreB(all domains)"] = frozenset(k for k in leaves if sub_of[k] == "MreB")
    groups["euk actin+Arps"] = EUK
    groups["euk actin+Arps+Asgard"] = EUK | ASG
    groups["crenactin+Ta0583+Asgard+euk"] = EUK | ASG | frozenset(
        k for k in leaves if sub_of[k] in ("crenactin", "Ta0583_archaeal_actin"))
else:
    groups["FtsZ all (bact+arch+euk organellar+Asgard)"] = frozenset(
        k for k in leaves if "FtsZ" in sub_of[k])
    groups["euk tubulins"] = EUK
    groups["euk tubulins+artubulin+BtubAB"] = EUK | frozenset(
        k for k in leaves if sub_of[k] in ("artubulin", "BtubA", "BtubB"))
    groups["euk_FtsZ+Cyanobacteria FtsZ"] = frozenset(
        k for k in leaves if sub_of[k] == "euk_FtsZ" or (sub_of[k] == "FtsZ" and sel.at[k, "group"] == "Cyanobacteriota"))
for g, mset in sorted(groups.items()):
    if len(mset) < 2:
        continue
    b = min_clade(mset)
    intr = sorted(b - mset)
    sup = np.nan
    for i, s in enumerate(df._side):
        if s == b or ALL - s == b:
            sup = df.ufboot.iloc[i]
    checks.append(dict(group=g, n=len(mset), is_split=len(intr) == 0, smallest_clade=len(b), n_intruders=len(intr),
                       ufboot=sup, intruders=";".join(f"{k}:{sub_of[k]}" for k in intr[:12])))
checks = pd.DataFrame(checks)
checks.to_csv(f"{OUT}/clade_checks_{sf}_{src}.tsv", sep="\t", index=False)
print(checks.drop(columns="intruders").to_string())

# minimal clade containing all eukaryotes (rooted on the prokaryotic side for this purpose)
mono = {"EUK_monophyletic": bool((df["class"] == "EUK_STEM").any()),
        "EUK+ASGARD_monophyletic": bool((df["class"] == "EUK+ASGARD_STEM").any()) if ASG else None}
# smallest split containing all EUK on one side: report the intruders
best = None
for s in df._side:
    for side in (s, ALL - s):
        if EUK <= side and (best is None or len(side) < len(best)):
            best = side
intr = sorted(best - EUK)
mono["smallest_clade_with_all_euk_size"] = len(best)
mono["non_euk_in_that_clade"] = [f"{k}:{dom_of[k]}:{sub_of[k]}" for k in intr][:40]
mono["n_non_euk_in_that_clade"] = len(intr)
mono["scale_sum3di_over_sumaa"] = float(scale)
mono["tree_length_aa"] = float(df.b_aa.sum()); mono["tree_length_3di"] = float(df.b_3di.sum())
fl = df[~df.terminal & df.supported & (df.b_aa >= 0.02)]
mono["n_supported_internal_baa_ge_0.02"] = int(len(fl))
mono["rho_supported_internal_quantiles"] = {q: float(fl.rho.quantile(q)) for q in (0.05, 0.25, 0.5, 0.75, 0.95)}
# focal splits (the 'eukaryote stem' and its nested alternatives)
SUBSET = lambda labs: frozenset(k for k in leaves if sub_of[k] in labs)
if sf == "actin":
    FOCAL = {"EUK_STEM": EUK, "EUK+ASGARD_STEM": EUK | ASG,
             "ASGARD_ACTINS_STEM": ASG,
             "EUK+ASGARD+crenactin+Ta0583_STEM": EUK | ASG | SUBSET({"crenactin", "Ta0583_archaeal_actin"})}
else:
    FOCAL = {"EUK_STEM": EUK, "TUBULIN_CLADE_STEM(euk+artubulin+BtubAB)": EUK | SUBSET({"artubulin", "BtubA", "BtubB"}),
             "TUBULIN+CetZ_STEM": EUK | SUBSET({"artubulin", "BtubA", "BtubB", "CetZ"}),
             "ARTUBULIN_STEM": SUBSET({"artubulin"}), "FtsZ_ALL_STEM": SUBSET({"FtsZ", "euk_FtsZ", "Asgard_FtsZ"})}
df["focal"] = ""
mono["focal"] = {}
for name, fs in FOCAL.items():
    if len(fs) < 1:
        continue
    hit = [i for i, sd in enumerate(df._side) if sd == fs or ALL - sd == fs]
    if not hit:
        mono["focal"][name] = "not a split in the ML tree"
        continue
    i = hit[0]
    df.loc[df.index[i], "focal"] = name
    r = df.iloc[i]
    mono["focal"][name] = dict(n_tips=len(fs), b_aa=float(r.b_aa), b_3di=float(r.b_3di), rho=float(r.rho),
                               ufboot=float(r.ufboot),
                               rho_percentile_among_supported_internal=float((fl.rho < r.rho).mean() * 100),
                               b_aa_percentile_among_supported_internal=float((fl.b_aa < r.b_aa).mean() * 100))
# eukaryote-internal supported branches vs prokaryote-internal supported branches
for c in ["within_eukaryote", "within_prok_or_mixed"]:
    x = fl[fl["class"] == c]
    mono[f"rho_{c}_supported_internal"] = dict(n=int(len(x)), median=float(x.rho.median()) if len(x) else None,
        q25=float(x.rho.quantile(.25)) if len(x) else None, q75=float(x.rho.quantile(.75)) if len(x) else None,
        pooled=float(x.b_3di.sum() / x.b_aa.sum() / scale) if len(x) else None)
st = df[df.clade_label != ""][["clade_label", "split_size", "ufboot", "b_aa", "b_3di", "rho"]]
mono["subfamily_stems"] = st.sort_values("clade_label").to_dict("records")
df = df.drop(columns="_side")
df.sort_values(["terminal", "class", "split_size"]).to_csv(f"{OUT}/branches_{sf}_{src}_{tag}.tsv", sep="\t",
                                                          index=False)
json.dump(mono, open(f"{OUT}/branch_summary_{sf}_{src}_{tag}.json", "w"), indent=1, default=str)
print(json.dumps({k: v for k, v in mono.items() if k != "subfamily_stems"}, indent=1, default=str))
print(st.sort_values("clade_label").to_string())
