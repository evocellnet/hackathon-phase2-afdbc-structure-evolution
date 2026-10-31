"""V11: independent E5 point quantities (§4.5 + D8): lineages, sigma2_obs, sigma2_L (= sigma2_obs - project floor F),
sister-contrast beta (c_v, vx_v from own HOG-pigeonhole replicates; several seeds)."""
import sys
import numpy as np
import dendropy
import scipy.sparse as sps
from v11_common import *

p4 = endpoint("E4"); p5 = endpoint("E5")
own4 = None
try:
    own4 = json.loads((OUT / "V11_branch.json").read_text())["E4"]["own"]
except Exception:
    pass
FE = {k: (own4 or p4)[k] for k in ("mu", "gamma", "x_center", "z_center")}
FE["beta"] = own4["beta"] if own4 else p4["estimate"]
TAU = TUN["tau_misfit"]
br = load_branches(PRIMARY, SPLIT)
fit = pl.read_parquet(ROOT / f"data/interim/branches/{PRIMARY}_fit.parquet")
f2 = fit.filter(pl.col("metric").is_in(["aa_ml", "tdi_ml"]))
bad = set(f2.filter(pl.col("excluded") | (pl.col("misfit_ratio").fill_nan(None) > TAU).fill_null(False))["hog"].to_list())
keep = sorted(set(fit["hog"].to_list()) - bad)
E = br.filter(pl.col("hog").is_in(keep))
A = E.filter(pl.col("metric") == "aa_ml").sort("hog", "edge_id"); T = E.filter(pl.col("metric") == "tdi_ml").sort("hog", "edge_id")
assert (A["edge_id"] == T["edge_id"]).all() and (A["hog"] == T["hog"]).all()
bA, bT = A["b"].to_numpy(), T["b"].to_numpy(); seT = T["b_se"].to_numpy()
hog = A["hog"].to_numpy(); child = A["panel_child"].to_numpy()
simple = A["panel_edge"].is_not_null().to_numpy()
rm = A["root_merged"].to_numpy()
res = {"n_hogs": len(keep), "root_merged_edges_with_panel_edge": int((simple & rm).sum())}

# ---------------- lineages
rbar = bT.sum() / bA.sum()
L = pl.DataFrame({"c": child[simple], "h": hog[simple], "a": bA[simple], "t": bT[simple], "v": seT[simple] ** 2,
                  "tau": A["tau"].to_numpy()[simple]}).group_by("c").agg(
    pl.col("h").n_unique().alias("nh"), pl.len().alias("n"), pl.col("a").sum().alias("D"), pl.col("t").sum().alias("S3"),
    pl.col("v").sum().alias("V"), pl.col("tau").min().alias("tau_min"), pl.col("tau").max().alias("tau_max")).sort("c")
L = L.with_columns((rbar * pl.col("D") / pl.col("V").sqrt()).alias("snr"))
Lk = L.filter((pl.col("nh") >= 50) & (pl.col("D") >= 0.05) & (pl.col("snr") >= TUN["e5"]["min_t_snr"]))
res["n_lineages_own"] = Lk.height; res["n_lineages_project"] = p5["n_lineages"]
res["n_lineages_without_D8"] = L.filter((pl.col("nh") >= 50) & (pl.col("D") >= 0.05)).height
assert np.allclose(Lk["tau_min"].to_numpy(), Lk["tau_max"].to_numpy())   # simple edges -> one panel edge -> one tau
S3 = Lk["S3"].to_numpy(); cens = S3 <= 0; S3 = np.where(cens, 0.005, S3)
D = Lk["D"].to_numpy()
logR = np.log(S3 / D)
x_e = np.log(Lk["tau_min"].to_numpy()) - FE["x_center"]
z_e = np.log(D / Lk["n"].to_numpy()) - FE["z_center"]
ell = logR - (FE["mu"] + FE["beta"] * x_e + FE["gamma"] * z_e)

# species pLDDT (pair tables, as the driver flag 'species_plddt_from_pair_tables'), descendant leaves from panel tree
pp = pairs().filter(pl.col("pair_set") == "primary")
sp = np.r_[pp["species_a"].to_numpy(), pp["species_b"].to_numpy()]; pv = np.r_[pp["mean_plddt_a"].to_numpy(), pp["mean_plddt_b"].to_numpy()]
spl = pl.DataFrame({"s": sp, "p": pv}).group_by("s").agg(pl.col("p").mean())
splddt = dict(zip(spl["s"].to_list(), spl["p"].to_list()))
tr = dendropy.Tree.get(path=str(ROOT / "data/interim/timetree/panel_tree.nwk"), schema="newick", preserve_underscores=True,
                       suppress_internal_node_taxa=True, suppress_leaf_node_taxa=False)
leaves_of, label_of, kids_of = {}, {}, {}
for nd in tr.postorder_node_iter():
    lab = str(int(nd.taxon.label.rsplit("_", 1)[-1])) if nd.is_leaf() else nd.label
    label_of[nd] = lab
    leaves_of[lab] = {int(lab)} if nd.is_leaf() else set().union(*[leaves_of[label_of[c]] for c in nd.child_nodes()])
    if not nd.is_leaf():
        kids_of[lab] = [label_of[c] for c in nd.child_nodes()]
pld = np.array([np.nanmean([splddt.get(t, np.nan) for t in leaves_of[str(c)]]) for c in Lk["c"].to_list()])
X = np.column_stack([np.ones(len(ell)), pld])
coef = np.linalg.lstsq(X, ell, rcond=None)[0]
rho_e = ell - X @ coef
s2obs = float(np.var(rho_e, ddof=1))
res["sigma2_obs_own"] = s2obs; res["sigma2_obs_project"] = p5["sigma2_obs"]; res["sigma2_obs_absdiff"] = abs(s2obs - p5["sigma2_obs"])
res["sigma2_L_own_with_project_F"] = s2obs - p5["floor_F"]; res["sigma2_L_project"] = p5["sigma2_L"]
res["sigma2_L_absdiff"] = abs(s2obs - p5["floor_F"] - p5["sigma2_L"])
res["n_censored_lineages"] = int(cens.sum())
print({k: v for k, v in res.items()}, flush=True)

# ---------------- sister contrasts
hogs_u, hcode = np.unique(hog, return_inverse=True)
H = len(hogs_u)
clades = [set(int(t) for t in e.split(",")) for e in A["edge_id"].to_list()]
contr = []
for v, ks in kids_of.items():
    ks = sorted(ks)
    for a, b in zip(ks[:-1], ks[1:]):
        contr.append((v, a, b))
starts = np.searchsorted(hcode, np.arange(H), "left"); ends = np.searchsorted(hcode, np.arange(H), "right")
rows, cols, ca, cb, meta = [], [], [], [], []
nrec = 0
for h in range(H):
    er = np.arange(starts[h], ends[h])
    cl = [clades[i] for i in er]
    S = set().union(*cl)
    root_merged = bool(rm[er].any())
    for ci, (v, a, b) in enumerate(contr):
        SA_ = S & leaves_of[a]; SB_ = S & leaves_of[b]
        if not SA_ or not SB_:
            continue
        need = SA_ | SB_
        cand = [j for j, c in enumerate(cl) if need <= c]
        if cand:
            u = cl[min(cand, key=lambda j: len(cl[j]))]
        else:
            if root_merged:
                continue
            u = S
        below = [j for j, c in enumerate(cl) if c <= u and len(c) < len(u)]
        for j in below:
            fa = len(cl[j] & SA_) / len(SA_); fb = len(cl[j] & SB_) / len(SB_)
            if fa or fb:
                rows.append(nrec); cols.append(er[j]); ca.append(fa); cb.append(fb)
        meta.append((ci, h)); nrec += 1
meta = np.array(meta)
CA = sps.csr_matrix((ca, (rows, cols)), shape=(nrec, len(bA))); CB = sps.csr_matrix((cb, (rows, cols)), shape=(nrec, len(bA)))
nh = np.bincount(meta[:, 0], minlength=len(contr))
kc = np.flatnonzero(nh >= 50); cmap = -np.ones(len(contr), int); cmap[kc] = np.arange(len(kc))
rk = cmap[meta[:, 0]]; m = rk >= 0
res["n_contrasts_own"] = len(kc); res["n_contrasts_project"] = p5["n_contrasts"]
LAa, LBa, LAt, LBt = (CA @ bA)[m], (CB @ bA)[m], (CA @ bT)[m], (CB @ bT)[m]
rh = meta[m, 1]; k = rk[m]


def xy(W):
    w = W[rh]
    with np.errstate(invalid="ignore", divide="ignore"):
        return (np.log(np.bincount(k, w * LAa, len(kc)) / np.bincount(k, w * LBa, len(kc))),
                np.log(np.bincount(k, w * LAt, len(kc)) / np.bincount(k, w * LBt, len(kc))))


x, y = xy(np.ones(H))
ok0 = np.isfinite(x) & np.isfinite(y)
res["n_contrasts_finite_xy_own"] = int(ok0.sum())
res["contrast_detail"] = [{"v": contr[c][0], "a": contr[c][1], "b": contr[c][2], "n_hogs": int(nh[c]), "x": float(x[i]), "y": float(y[i])}
                          for i, c in enumerate(kc)]
res["beta_naive_own"] = float(np.sum(x[ok0] * y[ok0]) / np.sum(x[ok0] ** 2))
betas = []
for seed in range(5):
    rng = np.random.default_rng([3331, seed])
    xs, ys = [], []
    for _ in range(200):
        W = rng.multinomial(H, np.full(H, 1 / H)).astype(float)
        a_, b_ = xy(W); xs.append(a_); ys.append(b_)
    xs, ys = np.array(xs), np.array(ys)
    c_v = np.nanmean((xs - np.nanmean(xs, 0)) * (ys - np.nanmean(ys, 0)), 0) * 200 / 199
    vx = np.nanvar(xs, 0, ddof=1)
    ok = ok0 & np.isfinite(c_v) & np.isfinite(vx)
    betas.append(float((np.sum(x[ok] * y[ok]) - np.sum(c_v[ok])) / (np.sum(x[ok] ** 2) - np.sum(vx[ok]))))
res["beta_own_seeds"] = betas; res["beta_own_mean"] = float(np.mean(betas)); res["beta_own_sd_over_seeds"] = float(np.std(betas, ddof=1))
res["beta_project"] = p5["beta_contrast"]; res["beta_absdiff"] = abs(np.mean(betas) - p5["beta_contrast"])
res["FE_source"] = "own E4 refit" if own4 else "project E4"
print(res, flush=True)
wjson(res, OUT / "V11_E5.json")
