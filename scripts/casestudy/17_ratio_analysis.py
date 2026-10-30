"""X6 node 2, steps 1-2: per-branch 3Di/AA ratio from IQ-TREE (fixed topology; 3 alignments; column bootstrap) and
from the project's branch engine (pairwise WLS on the UFBoot>=95 topology); flags and focal branches.

Branch identity = bipartition, oriented as the side NOT containing a fixed outgroup leaf (actin: an FtsA;
tubulin: a bacterial FtsZ) -> 'clade'. This orientation is for description only (no rooting claim).
rho_e = (b_3Di,e / b_AA,e) / (sum_e b_3Di / sum_e b_AA) per method (and per bootstrap replicate).
Comparable branches ('C'): supported (UFBoot >= 95 internal, or terminal), same kind (internal / terminal),
0.02 <= b_AA <= 2 (b_AA > 2 = 'saturated', excluded from flagging, reported).
Flag: log rho outside the 2.5-97.5 % range of the K = 50 C-branches nearest in log b_AA ('comparable length').
Bootstrap (columns resampled, same columns for AA and 3Di): rho 95 % percentile CI and p_beyond = fraction of
replicates beyond the same (full-data) window threshold. 'robust' = flagged with p_beyond >= 0.8 in the primary
(MAFFT L-INS-i) alignment and flagged in the same direction in >= 2 of the 3 alignments (linsi, famsa, structure-guided FoldMason).
Descriptive only (DESIGN 11.8). Usage: python 17_ratio_analysis.py <sf> [tree=ml|cons]
"""
import sys, os, glob, json
import numpy as np, pandas as pd
from ete3 import Tree

R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
CS = f"{R}/data/interim/casestudy"
OUT = f"{R}/results/casestudy"
sf = sys.argv[1]
TREE = sys.argv[2] if len(sys.argv) > 2 else "ml"
K, LO, HI = 50, 0.02, 2.0
sel = pd.read_csv(f"{OUT}/sample_members.tsv", sep="\t")
sel = sel[sel.superfamily == sf].set_index("seqid")
leaves = sorted(sel.index)
bit = {l: 1 << i for i, l in enumerate(leaves)}
ALL = (1 << len(leaves)) - 1
OG = sel[(sel.subfamily == ("FtsA" if sf == "actin" else "FtsZ")) & (sel.domain == "Bacteria")].index[0]


def orient(m):
    return (ALL ^ m) if (m & bit[OG]) else m


def members(m):
    return [l for l in leaves if m & bit[l]]


def tree_bl(path, support=False):
    t = Tree(path, format=1)
    out = {}
    for n in t.traverse("postorder"):
        if n.is_leaf():
            n.m = bit[n.name]
        else:
            n.m = 0
            for c in n.children:
                n.m |= c.m
        if n.is_root():
            continue
        m = orient(n.m)
        if m == 0 or m == ALL:
            continue
        sup = np.nan
        if support and not n.is_leaf():
            try:
                sup = float(str(n.name).split("/")[-1])
            except ValueError:
                pass
        if m in out:   # root-adjacent duplicate (rooted input): merge lengths
            out[m] = (out[m][0] + n.dist, max(out[m][1], sup) if support else sup)
        else:
            out[m] = (n.dist, sup)
    return out


# ---------------------------------------------------------------- topology & supports
ref_tree = f"{CS}/iqtree/{sf}_linsi_aa_ml.treefile" if TREE == "ml" else f"{CS}/constraint/{sf}_cons.treefile"
SUP = {m: s for m, (_, s) in tree_bl(f"{OUT}/{sf}_aa_ml_ufboot.nwk", support=True).items()}
TOP = tree_bl(ref_tree)
B = pd.DataFrame({"m": list(TOP)})
B["n_clade"] = [bin(m).count("1") for m in B.m]
B["terminal"] = (B.n_clade == 1) | (B.n_clade == len(leaves) - 1)
B["ufboot"] = [SUP.get(m, np.nan) for m in B.m]   # NaN for splits not in the ML tree (constrained tree only)
B["supported"] = B.terminal | (B.ufboot >= 95)


def describe(m):
    ms = members(m)
    s = sel.loc[ms]
    dom = s.domain.value_counts().to_dict()
    sub = s.subfamily.value_counts()
    grp = s.group.value_counts()
    d = {"domains": json.dumps(dom), "subfamilies": json.dumps(sub.head(5).to_dict()),
         "groups": json.dumps(grp.head(4).to_dict())}
    if len(ms) <= 3:
        d["taxa"] = "; ".join(f"{x}|{sel.at[x, 'name'].split(' (')[0]}|{sel.at[x, 'subfamily']}" for x in ms)
    else:
        d["taxa"] = ""
    return d


B = pd.concat([B, pd.DataFrame([describe(m) for m in B.m])], axis=1)

# ---------------------------------------------------------------- methods
SM = "fmg" if sf == "actin" else "foldmason"   # structure-guided MSA: actin = FoldMason guided by ML tree (default
# isolates FtsA/MreB rows); tubulin = default FoldMason (the guided one isolates rows instead)
srcs = ["linsi", "famsa", SM]
METH = {}
for src in srcs:
    pa = f"{CS}/iqfix/{sf}_{src}_{TREE}_aa_full.treefile"
    if not (os.path.exists(pa) and os.path.exists(pa.replace("_aa_full", "_3di_full"))):
        continue
    A = tree_bl(pa); T = tree_bl(f"{CS}/iqfix/{sf}_{src}_{TREE}_3di_full.treefile")
    bo = sorted(glob.glob(f"{CS}/iqfix/boot/{sf}_{src}_{TREE}_aa_b*.treefile"),
                key=lambda p: int(p.split("_b")[-1].split(".")[0]))
    reps = []
    for p in bo:
        q = p.replace("_aa_b", "_3di_b")
        if os.path.exists(q):
            a_, t_ = tree_bl(p), tree_bl(q)
            # IQ-TREE re-attaches sequences that are identical in a resampled alignment as cherries, which can move
            # a few near-zero splits (tubulin: RF 2-4 in ~1/3 of AA replicates); those splits are NaN in that replicate
            reps.append((np.array([a_.get(m, (np.nan,))[0] for m in B.m]), np.array([t_.get(m, (np.nan,))[0] for m in B.m])))
    METH[f"iq_{src}"] = (np.array([A[m][0] for m in B.m]), np.array([T[m][0] for m in B.m]), reps)
# LG+C20 AA lengths (saturation / LBA sensitivity of the AA side; 3Di side unchanged)
pc = f"{CS}/iqfix/{sf}_linsi_{TREE}_C20_aa_full.treefile"
if os.path.exists(pc):
    C = tree_bl(pc)
    B["b_aa_LG_C20"] = [C[m][0] for m in B.m]
# engine (collapsed topology): only bipartitions present in the reference topology
ENG = {}
if TREE == "ml":
    for v in ("tdiaa", "msa"):
        p = f"{OUT}/engine_branches_{sf}_{v}.tsv"
        if not os.path.exists(p):
            continue
        E = pd.read_csv(p, sep="\t")
        d = {}
        for _, r in E.iterrows():
            m = 0
            for x in r.leaves.split(";"):
                m |= bit[x]
            m = orient(m)
            if m in d:   # merged root edges: add shares back together
                d[m] = {k: d[m][k] + r[k] if k.startswith("b_") and not k.startswith("b_se") else r[k] for k in r.index}
            else:
                d[m] = r.to_dict()
        ENG[v] = d
        tot_aa = E.b_aa_ml.sum(); tot_td = E.b_tdi_ml.sum()
        B[f"eng_{v}_b_aa"] = [d.get(m, {}).get("b_aa_ml", np.nan) for m in B.m]
        B[f"eng_{v}_b_aa_se"] = [d.get(m, {}).get("b_se_aa_ml", np.nan) for m in B.m]
        B[f"eng_{v}_b_3di"] = [d.get(m, {}).get("b_tdi_ml", np.nan) for m in B.m]
        B[f"eng_{v}_b_3di_se"] = [d.get(m, {}).get("b_se_tdi_ml", np.nan) for m in B.m]
        cov = np.array([d.get(m, {}).get("cov_aa_tdi_ml", np.nan) for m in B.m])
        a_, t_ = B[f"eng_{v}_b_aa"].values, B[f"eng_{v}_b_3di"].values
        sa, st = B[f"eng_{v}_b_aa_se"].values, B[f"eng_{v}_b_3di_se"].values
        sc = tot_td / tot_aa
        with np.errstate(invalid="ignore", divide="ignore"):
            B[f"eng_{v}_rho"] = (t_ / a_) / sc
            # delta-method SE of log ratio (sum scale treated as fixed)
            B[f"eng_{v}_logrho_se"] = np.sqrt((st / t_) ** 2 + (sa / a_) ** 2 - 2 * cov / (t_ * a_))
        if v == "tdiaa":
            B["eng_tdiaa_b_lddtD"] = [d.get(m, {}).get("b_lddt_D", np.nan) for m in B.m]
            tl = E.b_lddt_D.sum()
            with np.errstate(invalid="ignore", divide="ignore"):
                B["eng_tdiaa_rho_lddt"] = (B.eng_tdiaa_b_lddtD / a_) / (tl / tot_aa)


def window_flags(b_aa, lrho, comp, kind):
    """2.5/97.5 % of log rho among the K comparable branches of the same kind nearest in log b_aa."""
    lo = np.full(len(b_aa), np.nan); hi = lo.copy(); med = lo.copy()
    lb = np.log(np.where(b_aa > 0, b_aa, np.nan))
    for i in range(len(b_aa)):
        if not np.isfinite(lb[i]) or not np.isfinite(lrho[i]):
            continue
        pool = np.flatnonzero(comp & (kind == kind[i]) & np.isfinite(lrho) & (np.arange(len(b_aa)) != i))
        if len(pool) < 20:
            continue
        nn = pool[np.argsort(np.abs(lb[pool] - lb[i]))[:K]]
        lo[i], med[i], hi[i] = np.quantile(lrho[nn], [0.025, 0.5, 0.975])
    return lo, med, hi


kind = B.terminal.values
for name, (a_, t_, reps) in METH.items():
    sc = t_.sum() / a_.sum()
    with np.errstate(invalid="ignore", divide="ignore"):
        rho = (t_ / a_) / sc
        lr = np.log(rho)
    comp = B.supported.values & (a_ >= LO) & (a_ <= HI)
    lo, med, hi = window_flags(a_, lr, comp, kind)
    B[f"{name}_b_aa"] = a_; B[f"{name}_b_3di"] = t_; B[f"{name}_rho"] = rho
    B[f"{name}_win_lo"] = np.exp(lo); B[f"{name}_win_med"] = np.exp(med); B[f"{name}_win_hi"] = np.exp(hi)
    flag = np.where(lr > hi, "high", np.where(lr < lo, "low", ""))
    flag = np.where(B.supported.values & (a_ >= LO) & (a_ <= HI), flag, np.where(a_ > HI, "saturated", ""))
    B[f"{name}_flag"] = flag
    if reps:
        RA = np.stack([r[0] for r in reps]); RT = np.stack([r[1] for r in reps])
        with np.errstate(invalid="ignore", divide="ignore"):
            RR = (RT / RA) / (np.nansum(RT, 1) / np.nansum(RA, 1))[:, None]
            LR = np.log(RR)
        B[f"{name}_rho_lo"] = np.nanpercentile(RR, 2.5, axis=0)
        B[f"{name}_rho_hi"] = np.nanpercentile(RR, 97.5, axis=0)
        B[f"{name}_b_aa_lo"] = np.nanpercentile(RA, 2.5, axis=0); B[f"{name}_b_aa_hi"] = np.nanpercentile(RA, 97.5, axis=0)
        B[f"{name}_p_beyond"] = np.where(flag == "high", np.nanmean(np.where(np.isfinite(LR), LR > hi[None], np.nan), 0),
                                         np.where(flag == "low", np.nanmean(np.where(np.isfinite(LR), LR < lo[None], np.nan), 0), np.nan))
        B[f"{name}_n_boot"] = len(reps)

if "iq_linsi" in METH:
    agree = np.zeros(len(B), int)
    for s in ("iq_linsi", "iq_famsa", f"iq_{SM}"):
        if f"{s}_flag" in B:
            agree += (B[f"{s}_flag"].values == B["iq_linsi_flag"].values) & np.isin(B[f"{s}_flag"].values, ["high", "low"])
    B["n_aln_agree"] = agree
    B["robust"] = np.isin(B.iq_linsi_flag, ["high", "low"]) & (B.get("iq_linsi_p_beyond", 1.0) >= 0.8) & (agree >= 2)

# ---------------------------------------------------------------- named clades (focal + subfamily best-F1 clades)
sub = sel.subfamily.to_dict(); dom = sel.domain.to_dict()
S = lambda f: sum(bit[l] for l in leaves if f(l))
EUK = S(lambda l: dom[l] == "Eukaryota" and sub[l] != "euk_FtsZ")
ASG = S(lambda l: dom[l] == "Asgard")
if sf == "actin":
    FOC = {"EUK stem (actin+Arps)": EUK, "EUK+Asgard stem": EUK | ASG, "Asgard actins": ASG,
           "crenactin+EUK+Asgard": EUK | ASG | S(lambda l: sub[l] == "crenactin"),
           "FtsA (all)": S(lambda l: sub[l] == "FtsA"), "MreB (bacterial)": S(lambda l: sub[l] == "MreB" and dom[l] == "Bacteria")}
else:
    FOC = {"EUK stem (euk tubulins)": EUK, "tubulin clade (euk+artubulin+BtubAB)": EUK | S(lambda l: sub[l] in ("artubulin", "BtubA", "BtubB")),
           "CetZ+tubulin": EUK | S(lambda l: sub[l] in ("artubulin", "BtubA", "BtubB", "CetZ")),
           "artubulins": S(lambda l: sub[l] == "artubulin"), "CetZ": S(lambda l: sub[l] == "CetZ")}
name_of = {}
for nm, m in FOC.items():
    for i, mm in enumerate(B.m):
        if mm == orient(m) or mm == orient(ALL ^ m):
            name_of[i] = nm
# subfamily clades: best F1 split per label
for lab in sorted(set(sub.values())):
    L = S(lambda l: sub[l] == lab)
    nL = bin(L).count("1")
    if nL < 3:
        continue
    best = (0, None)
    for i, mm in enumerate(B.m):
        if B.terminal.iat[i]:
            continue
        for side in (mm, ALL ^ mm):
            tp = bin(side & L).count("1")
            if tp == 0:
                continue
            f1 = 2 * tp / (bin(side).count("1") + nL)
            if f1 > best[0]:
                best = (f1, i)
    if best[0] >= 0.8 and best[1] not in name_of:
        name_of[best[1]] = f"{lab} stem (F1={best[0]:.2f})"
B["clade_name"] = [name_of.get(i, "") for i in range(len(B))]
B.insert(0, "branch_id", [f"{sf}_{TREE}_{i}" for i in range(len(B))])
B["clade_hex"] = [hex(m) for m in B.m]
B = B.drop(columns="m")
B.to_csv(f"{OUT}/branch_ratios_{sf}_{TREE}.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- summary
summ = {"tree": TREE, "n_branches": len(B), "outgroup_leaf_for_orientation": OG}
for name in METH:
    c = B.supported & (B[f"{name}_b_aa"] >= LO) & (B[f"{name}_b_aa"] <= HI)
    lr = np.log(B.loc[c, f"{name}_rho"])
    summ[name] = {"n_comparable_internal": int((c & ~B.terminal).sum()), "n_comparable_terminal": int((c & B.terminal).sum()),
                  "n_saturated_supported": int((B.supported & (B[f"{name}_b_aa"] > HI)).sum()),
                  "rho_q_internal": {q: float(np.exp(np.quantile(lr[~B.loc[c, "terminal"]], q))) for q in (0.025, 0.25, 0.5, 0.75, 0.975)},
                  "rho_q_terminal": {q: float(np.exp(np.quantile(lr[B.loc[c, "terminal"]], q))) for q in (0.025, 0.25, 0.5, 0.75, 0.975)},
                  "n_flag_high": int((B[f"{name}_flag"] == "high").sum()), "n_flag_low": int((B[f"{name}_flag"] == "low").sum()),
                  "tree_scale_3di_over_aa": float(B[f"{name}_b_3di"].sum() / B[f"{name}_b_aa"].sum())}
if "robust" in B:
    summ["n_robust_high"] = int((B.robust & (B.iq_linsi_flag == "high")).sum())
    summ["n_robust_low"] = int((B.robust & (B.iq_linsi_flag == "low")).sum())
cols = ["clade_name", "n_clade", "ufboot"] + [c for c in B.columns if c.endswith(("_b_aa", "_rho", "_rho_lo", "_rho_hi", "_flag", "_p_beyond"))]
foc = B[B.clade_name != ""]
summ["named"] = foc[cols].to_dict("records")
json.dump(summ, open(f"{OUT}/branch_ratio_summary_{sf}_{TREE}.json", "w"), indent=1, default=float)
pd.set_option("display.width", 250)
show = ["clade_name", "n_clade", "ufboot", "iq_linsi_b_aa", "iq_linsi_rho", "iq_linsi_rho_lo", "iq_linsi_rho_hi",
        "iq_linsi_flag", "iq_famsa_rho", "iq_fmg_rho", "eng_tdiaa_b_aa", "eng_tdiaa_rho", "eng_msa_rho"]
print(foc[[c for c in show if c in foc]].round(3).to_string())
print(json.dumps({k: v for k, v in summ.items() if k != "named"}, indent=1, default=float))

# ---------------------------------------------------------------- stem paths (adjacent short branches trade length)
# rho of the summed path of consecutive named stem branches, full data and per bootstrap replicate
PATHS = {"actin": {"path: below crenactin+EUK+Asgard -> EUK+Asgard ancestor": ["crenactin+EUK+Asgard", "EUK+Asgard stem"]},
         "tubulin": {"path: CetZ+tubulin -> tubulin clade ancestor": ["CetZ+tubulin", "tubulin clade (euk+artubulin+BtubAB)"]}}[sf]
if TREE == "cons":
    PATHS = {"actin": {"path: below crenactin+EUK+Asgard -> EUK ancestor": ["crenactin+EUK+Asgard", "EUK+Asgard stem", "EUK stem (actin+Arps)"]},
             "tubulin": {"path: CetZ+tubulin -> EUK tubulin ancestor": ["CetZ+tubulin", "tubulin clade (euk+artubulin+BtubAB)", "EUK stem (euk tubulins)"]}}[sf]
prow = []
for pname, parts in PATHS.items():
    idx = [int(np.flatnonzero(B.clade_name.values == p)[0]) for p in parts if (B.clade_name.values == p).any()]
    for name, (a_, t_, reps) in METH.items():
        sc = t_.sum() / a_.sum()
        r = dict(sf=sf, tree=TREE, path=pname, parts=";".join(B.clade_name.values[i] for i in idx), method=name,
                 b_aa=float(a_[idx].sum()), b_3di=float(t_[idx].sum()), rho=float(t_[idx].sum() / a_[idx].sum() / sc))
        if reps:
            rr = np.array([(rt[idx].sum() / ra[idx].sum()) / (np.nansum(rt) / np.nansum(ra)) for ra, rt in reps])
            r["rho_lo"], r["rho_hi"] = np.percentile(rr, [2.5, 97.5])
            # percentile of the path rho among comparable internal branches of similar length (single branches)
        comp = B.supported.values & ~B.terminal.values & (a_ >= LO) & (a_ <= HI)
        lb = np.log(a_[comp]); lr_ = np.log((t_[comp] / a_[comp]) / sc)
        nn = np.argsort(np.abs(lb - np.log(r["b_aa"])))[:K]
        r["pct_in_window"] = float((lr_[nn] < np.log(r["rho"])).mean() * 100)
        r["pct_all_internal"] = float((lr_ < np.log(r["rho"])).mean() * 100)
        prow.append(r)
    for v in ("tdiaa", "msa"):
        if f"eng_{v}_b_aa" in B:
            a_ = B[f"eng_{v}_b_aa"].values; t_ = B[f"eng_{v}_b_3di"].values
            if np.all(np.isfinite(a_[idx])):
                E = pd.read_csv(f"{OUT}/engine_branches_{sf}_{v}.tsv", sep="\t")
                sc = E.b_tdi_ml.sum() / E.b_aa_ml.sum()
                prow.append(dict(sf=sf, tree=TREE, path=pname, parts=";".join(B.clade_name.values[i] for i in idx),
                                 method=f"eng_{v}", b_aa=float(a_[idx].sum()), b_3di=float(t_[idx].sum()),
                                 rho=float(t_[idx].sum() / a_[idx].sum() / sc)))
if "b_aa_LG_C20" in B and "iq_linsi" in METH:
    # LG+C20+G4 AA lengths with the LG+G4-tree 3Di lengths (3Di model unchanged): rho under a site-heterogeneous AA model
    a20 = B.b_aa_LG_C20.values; t_ = METH["iq_linsi"][1]
    sc20 = t_.sum() / a20.sum()
    with np.errstate(invalid="ignore", divide="ignore"):
        B["iq_linsi_rho_C20"] = (t_ / a20) / sc20
    for pname, parts in PATHS.items():
        idx = [int(np.flatnonzero(B.clade_name.values == p)[0]) for p in parts if (B.clade_name.values == p).any()]
        prow.append(dict(sf=sf, tree=TREE, path=pname, parts=";".join(B.clade_name.values[i] for i in idx),
                         method="iq_linsi_LG+C20(AA)", b_aa=float(a20[idx].sum()), b_3di=float(t_[idx].sum()),
                         rho=float(t_[idx].sum() / a20[idx].sum() / sc20)))
    B.to_csv(f"{OUT}/branch_ratios_{sf}_{TREE}.tsv", sep="\t", index=False)
P = pd.DataFrame(prow)
P.to_csv(f"{OUT}/stem_paths_{sf}_{TREE}.tsv", sep="\t", index=False)
print(P.round(3).to_string())
