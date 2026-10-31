"""V11: independent D31 parametric-bootstrap bias of SIMEX for E1/E2 (approximate: own RNG stream, K reps).
Truth design = tree-smoothed AA distance x_hat (sum of AA b_raw along the HOG-tree path, from the branch table);
pseudo data: x* = x_hat + sx z1, y* = eta + a_h(c_S) x_hat^c_S + sy(rho z1 + sqrt(1-rho^2) z2) + sqrt(s2) z3;
re-estimated by two-pass + SIMEX (B_s=20). Usage: python v11_curve_pb.py E1 K c_simex_ref"""
import sys, json, time
import numpy as np
from multiprocessing import Pool
from v11_common import *
import v11_curve as vc

EP = sys.argv[1]; K = int(sys.argv[2])
vc.EP = EP; vc.CFG = {"E1": dict(y="tdi_ml", se="tdi_ml_se_k", sat="tdi_sat", eta=TUN["eta"][PRIMARY]["tdi_ml"], rho="aa_tdi", ff="tdi"),
                      "E2": dict(y="lddt_D", se="lddt_D_se", sat="lddt_sat", eta=TUN["eta"][PRIMARY]["lddt_D"], rho="aa_lddt", ff="lddt")}[EP]
CFG = vc.CFG
d, cnt, _ = vc.prepare()
hog = d["hog"].to_numpy(); x = d["aa_ml"].to_numpy().astype(float); y = d[CFG["y"]].to_numpy().astype(float)
sy = d[CFG["se"]].to_numpy().astype(float); sx = d["aa_ml_se_k"].to_numpy().astype(float)
rt = TUN["rho"][PRIMARY]; rho = binlook(rt["bin_edges"], rt[CFG["rho"]], x)
F = vc.Fit(hog, x, y, CFG["eta"])

# ---- tree-smoothed x_hat from the AA branch table
br = load_branches(PRIMARY, SPLIT).filter(pl.col("metric") == "aa_ml")
fit = pl.read_parquet(ROOT / f"data/interim/branches/{PRIMARY}_fit.parquet").filter(pl.col("metric") == "aa_ml")
under = set(fit.filter(pl.col("reason").str.contains("underdetermined|rank_deficient"))["hog"].to_list())
edges = {}
for h, e, b in br.select("hog", "edge_id", "b_raw").iter_rows():
    edges.setdefault(h, []).append((set(int(t) for t in e.split(",")), b))
sa = d["species_a"].to_numpy(); sb = d["species_b"].to_numpy()
xhat = np.full(len(x), np.nan)
for i in range(len(x)):
    h = hog[i]
    if h in under or h not in edges:
        continue
    s = 0.0
    for clade, b in edges[h]:
        if (sa[i] in clade) != (sb[i] in clade):
            s += b
    xhat[i] = s
design_ok = np.isfinite(xhat) & (xhat > 0)
xd = np.where(design_ok, xhat, x)
print("x_hat: n_nofit/nonpos", int((~design_ok).sum()), "median |xhat-x|/x", float(np.median(np.abs(xhat[design_ok] - x[design_ok]) / x[design_ok])), flush=True)

c1, _ = F.chat(1 / sy ** 2); _, _, r1 = F.prof(c1, 1 / sy ** 2)
s2 = max(0.0, np.mean(r1 ** 2) - np.mean(sy ** 2)); w = 1 / (sy ** 2 + s2)
c0, _ = F.chat(w)
cS = float(sys.argv[3])     # SIMEX estimate (own run) used as pseudo-truth
_, a, _ = F.prof(cS, w)
mu_c = a[F.g] * np.exp(cS * np.log(xd))     # y - eta mean
LAM = [0.5, 1.0, 1.5, 2.0]


def one(k):
    rng = np.random.default_rng([2211, k])
    z1, z2, z3 = rng.standard_normal((3, len(x)))
    xs = np.maximum(xd + sx * z1, 1e-4)
    ys = mu_c + sy * (rho * z1 + np.sqrt(1 - rho ** 2) * z2) + np.sqrt(s2) * z3
    Fk = vc.Fit(hog, xs, ys + CFG["eta"], CFG["eta"])
    win = (cS - 0.3, cS + 0.3)
    ck1, _ = Fk.chat(1 / sy ** 2, window=win); _, _, rk = Fk.prof(ck1, 1 / sy ** 2)
    s2k = max(0.0, np.mean(rk ** 2) - np.mean(sy ** 2)); wk = 1 / (sy ** 2 + s2k)
    ck, _ = Fk.chat(wk, window=win)
    cl = {0.0: ck}
    for lam in LAM:
        cs = []
        for b in range(20):
            u1 = rng.standard_normal(len(x)); u2 = rng.standard_normal(len(x))
            xb = np.maximum(xs + np.sqrt(lam) * sx * u1, 1e-4)
            yb = Fk.yc + np.sqrt(lam) * sy * (rho * u1 + np.sqrt(1 - rho ** 2) * u2)
            cs.append(Fk.chat(wk, yc=yb, lx=np.log(xb), window=(ck - 0.3, ck + 0.3))[0])
        cl[lam] = float(np.mean(cs))
    L = np.array(sorted(cl)); C = np.array([cl[l] for l in L])
    return float(np.polyval(np.polyfit(L, C, 2), -1.0))


if __name__ == "__main__":
    t0 = time.time()
    with Pool(min(K, 25)) as p:
        reps = np.array(p.map(one, range(K)))
    bias = float(reps.mean() - cS); mc = float(reps.std(ddof=1) / np.sqrt(K))
    proj = endpoint(EP)
    out = {"endpoint": EP, "K": K, "c_S_used": cS, "reps": reps.tolist(), "pb_bias_own": bias, "pb_mc_se_own": mc,
           "pb_bias_project": proj["pb_bias"], "pb_mc_se_project": proj["pb_mc_se"],
           "pb_bias_absdiff": abs(bias - proj["pb_bias"]), "z_diff": abs(bias - proj["pb_bias"]) / np.hypot(mc, proj["pb_mc_se"]),
           "n_design_nofit_or_nonpos": int((~design_ok).sum()), "runtime_s": time.time() - t0}
    print(out)
    wjson(out, OUT / f"V11_{EP}_pb.json")
