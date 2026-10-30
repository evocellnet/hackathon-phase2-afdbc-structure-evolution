"""X7 step 8: per-branch b_3Di / b_AA on the fixed AA ML topology (descriptive).

rho_e = (b_3Di,e / b_AA,e) / (sum b_3Di / sum b_AA) over all branches of that run (3Di and AA scales differ).
Rooting (for classifying clades only; branch lengths are unrooted): on the clade = MRCA of all ELK/atypical domains
(RIO1, ABC1, APH, Choline_kinase, Kdo), computed with the tree first rooted at a eukaryotic ePK tip; intruders listed.
Coarse categories of tips: euk_ePK (eukaryotic Pkinase/PK_Tyr), prok_ePK (bacterial/archaeal/Asgard Pkinase/PK_Tyr),
euk_ELK, prok_ELK. Branch classes: terminal_<cat>; internal 'within_<cat>' if the clade below is pure, else 'mixed'.
'supported' = terminal or UFBoot >= 95 (from the ML tree). Background sets use supported internal branches with
b_AA >= 0.02 (as X6). Group stems: for each eukaryotic ePK group, the MRCA clade of its members (purity reported)
and the largest pure clade (>= 3 tips); for TK also the MRCA of holozoan (Metazoa + Holozoa_other) TKs.
Flags: 'rho_high' = supported internal with rho > q95 of the supported-internal background; 'rho_low' < q05;
'sat_risk' = b_AA > 1.0 subs/site (LBA/saturation-prone).
Usage: python 08_branch_ratio.py <src>
"""
import sys, json
import numpy as np, pandas as pd
from ete3 import Tree

R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
K = f"{R}/data/interim/kinase"
OUT = f"{R}/results/kinase"
src = sys.argv[1]
ml = Tree(f"{K}/iqtree/pkl_linsi_aa_ml.treefile", format=1)
aa = Tree(f"{K}/iqtree/pkl_{src}_aa_fixed.treefile", format=1)
td = Tree(f"{K}/iqtree/pkl_{src}_3di_fixed.treefile", format=1)
sel = pd.read_csv(f"{OUT}/sample_members.tsv", sep="\t").set_index("seqid")
ELKF = ("RIO1", "ABC1", "APH", "Choline_kinase", "Kdo")


def cat(k):
    c = sel.at[k, "klass"]
    if c.startswith("ePK_"):
        return "euk_ePK"
    if c in ("prok_ePK", "Asgard_ePK"):
        return "prok_ePK"
    return "euk_ELK" if c.endswith("_euk") else "prok_ELK"


CAT = {k: cat(k) for k in sel.index}
GRP = {k: (sel.at[k, "klass"][4:] if CAT[k] == "euk_ePK" else "") for k in sel.index}
HOLO = {k for k in sel.index if sel.at[k, "group"] in ("Metazoa", "Holozoa_other")}
leaves = ml.get_leaf_names()
ALL = frozenset(leaves)


def sup_of(n):
    try:
        return float(str(n.name).split("/")[-1])
    except ValueError:
        return np.nan


# --- root on the ELK clade
anyepk = next(k for k in leaves if CAT[k] == "euk_ePK")
ml.set_outgroup(ml & anyepk)
elk = [k for k in leaves if CAT[k] in ("euk_ELK", "prok_ELK")]
mr = ml.get_common_ancestor(elk)
elk_intr = [k for k in mr.get_leaf_names() if k not in set(elk)]
ml.set_outgroup(mr)


def bip(s):
    s = frozenset(s)
    return s if leaves[0] not in s else ALL - s


Laa = {bip(n.get_leaf_names()): n.dist for n in aa.traverse() if not n.is_root()}
Ltd = {bip(n.get_leaf_names()): n.dist for n in td.traverse() if not n.is_root()}
rows, seen = [], set()
for n in ml.traverse("preorder"):
    if n.is_root():
        continue
    cl = n.get_leaf_names()
    b = bip(cl)
    if b in seen or len(b) in (0, len(ALL)):
        continue
    seen.add(b)
    cats = pd.Series([CAT[k] for k in cl]).value_counts()
    term = n.is_leaf()
    klass = (f"terminal_{cats.index[0]}" if term else (f"within_{cats.index[0]}" if len(cats) == 1 else "mixed"))
    grps = pd.Series([GRP[k] for k in cl if GRP[k]]).value_counts()
    rows.append(dict(node=n, clade_size=len(cl), terminal=term, ufboot=np.nan if term else sup_of(n),
                     b_aa=Laa.get(b, np.nan), b_3di=Ltd.get(b, np.nan), bclass=klass,
                     tip=cl[0] if term else "", cats=json.dumps(cats.to_dict()),
                     egroups=json.dumps(grps.head(5).to_dict()),
                     frac_holozoa=np.mean([k in HOLO for k in cl])))
df = pd.DataFrame(rows)
assert df.b_aa.notna().all() and df.b_3di.notna().all(), "topology mismatch"
# root edge: the two root children are one unrooted branch; the bip() dedup keeps the first (the ELK side)
scale = df.b_3di.sum() / df.b_aa.sum()
df["rho"] = df.b_3di / df.b_aa / scale
df["log2rho"] = np.log2(df.rho.clip(lower=1e-3))
df["supported"] = df.terminal | (df.ufboot >= 95)
bg = df[~df.terminal & df.supported & (df.b_aa >= 0.02)]
q05, q95 = bg.rho.quantile(.05), bg.rho.quantile(.95)
df["flag"] = ""
m_int = ~df.terminal & df.supported & (df.b_aa >= 0.02)
df.loc[m_int & (df.rho > q95), "flag"] += "rho_high;"
df.loc[m_int & (df.rho < q05), "flag"] += "rho_low;"
df.loc[df.b_aa > 1.0, "flag"] += "sat_risk;"
pct = lambda v, ref: float((ref < v).mean() * 100)

summ = dict(src=src, n_tips=len(ALL), scale_sum3di_over_sumaa=float(scale), tree_len_aa=float(df.b_aa.sum()),
            tree_len_3di=float(df.b_3di.sum()), elk_mrca_intruders=[f"{k}:{sel.at[k, 'klass']}" for k in elk_intr],
            n_supported_internal_bg=int(len(bg)), rho_bg_q05=float(q05), rho_bg_q95=float(q95),
            n_internal=int((~df.terminal).sum()), n_internal_ufb95=int((~df.terminal & df.supported).sum()))
# --- class distributions (question ii)
cls = []
for c, g in df.groupby("bclass"):
    for lab, x in (("all", g), ("supported_baa>=0.02", g[g.supported & (g.b_aa >= 0.02)])):
        if len(x) == 0:
            continue
        cls.append(dict(bclass=c, subset=lab, n=len(x), median_rho=x.rho.median(), q25=x.rho.quantile(.25),
                        q75=x.rho.quantile(.75), pooled_rho=x.b_3di.sum() / x.b_aa.sum() / scale,
                        median_b_aa=x.b_aa.median()))
cls = pd.DataFrame(cls)
cls.to_csv(f"{OUT}/class_rho_{src}.tsv", sep="\t", index=False)
# b_AA-matched comparison: median rho per class within b_AA quantile bins (supported branches, b_AA >= 0.02)
s = df[df.supported & (df.b_aa >= 0.02)].copy()
s["baa_bin"] = pd.qcut(s.b_aa, 5)
s["kind"] = np.where(s.terminal, "terminal", "internal")
s["cat"] = s.bclass.str.replace("terminal_", "").str.replace("within_", "")
binned = s.groupby(["kind", "baa_bin", "cat"], observed=True).rho.agg(["size", "median"]).reset_index()
binned["baa_bin"] = binned.baa_bin.astype(str)
binned.to_csv(f"{OUT}/class_rho_baa_binned_{src}.tsv", sep="\t", index=False)
# --- group stems (question i)
bg_euk = df[(df.bclass == "within_euk_ePK") & df.supported & (df.b_aa >= 0.02)]
stems = []
node_row = {id(r.node): i for i, r in df.iterrows()}
bip_row = {bip(r.node.get_leaf_names()): i for i, r in df.iterrows()}


def stem(name, members):
    members = [k for k in members if k in ALL]
    if len(members) < 2:
        return
    mrca = ml.get_common_ancestor(members)
    if mrca.is_root():
        stems.append(dict(stem=name, n_members=len(members), note="MRCA is root")); return
    i = node_row.get(id(mrca))
    note = ""
    if i is None:  # second root child: same unrooted branch as the first root child
        i = bip_row[bip(mrca.get_leaf_names())]
        note = "root-adjacent (clade = everything outside the ELK-MRCA clade)"
    r = df.loc[i]
    csize = len(mrca.get_leaf_names())
    stems.append(dict(stem=name, n_members=len(members), clade_size=csize,
                      purity=len(members) / csize, ufboot=r.ufboot, b_aa=r.b_aa, b_3di=r.b_3di, rho=r.rho,
                      rho_pct_vs_within_euk_ePK=pct(r.rho, bg_euk.rho), rho_pct_vs_all_bg=pct(r.rho, bg.rho),
                      b_aa_pct_vs_all_bg=pct(r.b_aa, bg.b_aa), flag=r.flag, note=note))


def largest_pure(name, members):
    ms = set(members)
    best = None
    for n in ml.traverse():
        if n.is_root() or n.is_leaf():
            continue
        lv = n.get_leaf_names()
        if len(lv) >= 3 and set(lv) <= ms and (best is None or len(lv) > len(best.get_leaf_names())):
            best = n
    if best is None:
        return
    r = df.loc[node_row[id(best)]]
    stems.append(dict(stem=name + " (largest pure clade)", n_members=len(ms), clade_size=r.clade_size,
                      purity=1.0, ufboot=r.ufboot, b_aa=r.b_aa, b_3di=r.b_3di, rho=r.rho,
                      rho_pct_vs_within_euk_ePK=pct(r.rho, bg_euk.rho), rho_pct_vs_all_bg=pct(r.rho, bg.rho),
                      b_aa_pct_vs_all_bg=pct(r.b_aa, bg.b_aa), flag=r.flag, note=""))


for g in ["AGC", "CAMK", "CK1", "CMGC", "STE", "TK", "TKL", "Other"]:
    mem = [k for k in leaves if GRP[k] == g]
    stem(f"{g} MRCA", mem); largest_pure(g, mem)
stem("TK holozoan MRCA", [k for k in leaves if GRP[k] == "TK" and k in HOLO])
largest_pure("TK holozoan", [k for k in leaves if GRP[k] == "TK" and k in HOLO])
stem("euk ePK MRCA", [k for k in leaves if CAT[k] == "euk_ePK"])
stem("prok ePK MRCA", [k for k in leaves if CAT[k] == "prok_ePK"])
largest_pure("prok ePK", [k for k in leaves if CAT[k] == "prok_ePK"])
stem("all ePK (euk+prok) MRCA", [k for k in leaves if CAT[k] in ("euk_ePK", "prok_ePK")])
for f in ELKF:
    mem = [k for k in leaves if sel.at[k, "klass"].startswith(f + "_")]
    stem(f"{f} MRCA", mem); largest_pure(f, mem)
stems = pd.DataFrame(stems)
stems.to_csv(f"{OUT}/group_stems_{src}.tsv", sep="\t", index=False)
summ["within_euk_ePK_bg"] = dict(n=int(len(bg_euk)), median=float(bg_euk.rho.median()),
                                 q25=float(bg_euk.rho.quantile(.25)), q75=float(bg_euk.rho.quantile(.75)))
# group purity / monophyly
mono = []
for g in ["AGC", "CAMK", "CK1", "CMGC", "STE", "TK", "TKL", "Other"]:
    mem = [k for k in leaves if GRP[k] == g]
    mono.append(dict(group=g, n=len(mem)))
summ["flagged_rho_high"] = int(df.flag.str.contains("rho_high").sum())
summ["flagged_sat_risk_internal"] = int((~df.terminal & df.flag.str.contains("sat_risk")).sum())
out = df.drop(columns="node")
out.to_csv(f"{OUT}/branches_{src}.tsv", sep="\t", index=False)
json.dump(summ, open(f"{OUT}/branch_summary_{src}.json", "w"), indent=1, default=float)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
print(json.dumps(summ, indent=1, default=float))
print(cls.to_string())
print(stems.to_string())
print(out[out.flag.str.contains("rho_high")][["clade_size", "ufboot", "b_aa", "b_3di", "rho", "bclass", "cats",
                                                "egroups", "flag"]].to_string())
