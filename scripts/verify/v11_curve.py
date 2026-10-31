"""V11: independent E1/E2 point estimates (naive c-hat, SIMEX, D31 parametric-bootstrap bias) from §4.1/§4.2.
Usage: python v11_curve.py [E1|E2] [--simex] [--pb K]"""
import sys, time
import numpy as np
from scipy.optimize import minimize_scalar
from v11_common import *

EP = sys.argv[1] if len(sys.argv) > 1 else "E1"
CFG = {"E1": dict(y="tdi_ml", se="tdi_ml_se_k", sat="tdi_sat", eta=TUN["eta"][PRIMARY]["tdi_ml"], rho="aa_tdi", ff="tdi"),
       "E2": dict(y="lddt_D", se="lddt_D_se", sat="lddt_sat", eta=TUN["eta"][PRIMARY]["lddt_D"], rho="aa_lddt", ff="lddt")}[EP]


def fail_frac_check(df):
    """Own fail_frac per HOG (primary pairs; not ok or saturated on the metric or AA) vs hogs_endpoint_cols."""
    st = df["status"].to_numpy() != "ok"
    sat = df[CFG["sat"]].fill_null(False).to_numpy().astype(bool) | df["aa_sat"].fill_null(False).to_numpy().astype(bool)
    y = df[CFG["y"]].to_numpy().astype(float)
    bad = st | sat | ~np.isfinite(y)
    t = pl.DataFrame({"hog": df["hog"], "bad": bad}).group_by("hog").agg(pl.col("bad").mean().alias("ff"))
    hc = hogcols().select("hog", f"fail_frac_{CFG['ff']}_{PRIMARY}", f"excl_{CFG['ff']}_{PRIMARY}")
    m = t.join(hc, on="hog", how="left")
    d = np.abs(m["ff"].to_numpy() - m[f"fail_frac_{CFG['ff']}_{PRIMARY}"].to_numpy())
    return {"n_hogs": m.height, "max_abs_diff_fail_frac": float(np.nanmax(d)),
            "n_excl_own": int((m["ff"] > 0.10).sum()), "n_excl_project": int(m[f"excl_{CFG['ff']}_{PRIMARY}"].sum())}, \
        set(t.filter(pl.col("ff") > 0.10)["hog"].to_list())


def prepare():
    df = pairs().filter(pl.col("pair_set") == "primary")
    n0 = df.height
    ffc, excl = fail_frac_check(df)
    cnt = {"primary_pairs": n0}
    one = df.filter((pl.col("one_to_one") == 1) & (pl.col("relation") == "ortholog"))
    cnt["not_1to1"] = n0 - one.height
    r = one.filter((pl.col("aa_ml_seq") >= 0.05) & (pl.col("aa_ml_seq") <= 2.0))
    cnt["aa_ml_seq_out_of_range"] = one.height - r.height
    r2 = r.filter(~pl.col("hog").is_in(list(excl)))
    cnt["hog_excluded_metric"] = r.height - r2.height
    x = r2["aa_ml"].to_numpy().astype(float); y = r2[CFG["y"]].to_numpy().astype(float)
    se = r2[CFG["se"]].to_numpy().astype(float); sx = r2["aa_ml_se_k"].to_numpy().astype(float)
    sat = r2[CFG["sat"]].fill_null(False).to_numpy().astype(bool) | r2["aa_sat"].fill_null(False).to_numpy().astype(bool)
    elig = (r2["status"].to_numpy() == "ok") & ~sat & np.isfinite(x) & np.isfinite(y) & np.isfinite(se) & (se > 0) & (x > 0)
    cnt["pair_weight0_fail_or_sat"] = int((~elig).sum())
    r3 = r2.filter(pl.Series(elig))
    g = r3.group_by("hog").agg(pl.len().alias("n"), (pl.col("aa_ml").max() / pl.col("aa_ml").min()).alias("rng"))
    keep = set(g.filter((pl.col("n") >= 10) & (pl.col("rng") >= 4))["hog"].to_list())
    r4 = r3.filter(pl.col("hog").is_in(list(keep)))
    cnt["hog_lt_min_pairs_or_range"] = r3.height - r4.height
    cnt["n_fit"] = r4.height; cnt["n_hogs"] = len(keep)
    return r4, cnt, ffc


class Fit:
    def __init__(self, hog, x, y, eta):
        _, self.g = np.unique(hog, return_inverse=True)
        self.H = self.g.max() + 1
        self.lx = np.log(x); self.yc = y - eta

    def prof(self, c, w, yc=None, lx=None):
        yc = self.yc if yc is None else yc; lx = self.lx if lx is None else lx
        xc = np.exp(c * lx)
        num = np.bincount(self.g, w * yc * xc, self.H); den = np.bincount(self.g, w * xc * xc, self.H)
        a = num / den
        r = yc - a[self.g] * xc
        return float(np.sum(w * r * r)), a, r

    def chat(self, w, yc=None, lx=None, window=None):
        grid = np.round(np.arange(0.2, 3.0 + 1e-9, 0.01), 2)
        if window is not None:   # speed-up for Monte-Carlo refits: search only c0 +- 0.3 (boundary hit is reported)
            grid = grid[(grid >= window[0] - 1e-9) & (grid <= window[1] + 1e-9)]
        Q = np.array([self.prof(c, w, yc, lx)[0] for c in grid])
        i = int(np.argmin(Q))
        lo, hi = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
        o = minimize_scalar(lambda c: self.prof(c, w, yc, lx)[0], bounds=(lo, hi), method="bounded",
                            options={"xatol": 1e-10})
        c = o.x if o.fun <= Q[i] else grid[i]
        return float(c), bool(i in (0, len(grid) - 1))


def main():
    t0 = time.time()
    d, cnt, ffc = prepare()
    proj = endpoint(EP)
    hog = d["hog"].to_numpy(); x = d["aa_ml"].to_numpy().astype(float); y = d[CFG["y"]].to_numpy().astype(float)
    sy = d[CFG["se"]].to_numpy().astype(float); sx = d["aa_ml_se_k"].to_numpy().astype(float)
    F = Fit(hog, x, y, CFG["eta"])
    c1, _ = F.chat(1 / sy ** 2)
    _, _, r1 = F.prof(c1, 1 / sy ** 2)
    s2 = max(0.0, np.mean(r1 ** 2) - np.mean(sy ** 2))
    w = 1 / (sy ** 2 + s2)
    c0, bnd = F.chat(w)
    out = {"endpoint": EP, "counts_own": cnt, "counts_project": proj["n_excluded"] | {"n_units": proj["n_units"], "n_hogs": proj["n_hogs"]},
           "fail_frac_check": ffc, "s2_own": s2, "s2_project": proj["s2"], "c_naive_own": c0, "c_naive_project": proj["c_naive"],
           "c_naive_absdiff": abs(c0 - proj["c_naive"]), "grid_boundary": bnd}
    print(out, f"{time.time()-t0:.0f}s", flush=True)
    if "--simex" in sys.argv:
        # §4.1 SIMEX with per-pair rho(bin) (D31 text: "now with per-pair rho(bin)"); kappa already in SEs.
        rt = TUN["rho"][PRIMARY]
        rho = binlook(rt["bin_edges"], rt[CFG["rho"]], x)
        rng = np.random.default_rng(1110)  # own seed (not the project's stream)
        Bs = 50; lams = [0.5, 1.0, 1.5, 2.0]; cl = {0.0: c0}; nclip = 0
        for lam in lams:
            cs = []
            for b in range(Bs):
                z1 = rng.standard_normal(len(x)); z2 = rng.standard_normal(len(x))
                ex = np.sqrt(lam) * sx * z1
                ey = np.sqrt(lam) * sy * (rho * z1 + np.sqrt(1 - rho ** 2) * z2)
                xb = x + ex; nclip += int((xb < 1e-4).sum()); xb = np.maximum(xb, 1e-4)
                cs.append(F.chat(w, yc=F.yc + ey, lx=np.log(xb))[0])   # pass-2 weights held (D5)
            cl[lam] = float(np.mean(cs)); print(lam, cl[lam], np.std(cs) / np.sqrt(Bs), flush=True)
        L = np.array(sorted(cl)); C = np.array([cl[l] for l in L])
        q = np.polyfit(L, C, 2); c_simex = float(np.polyval(q, -1.0))
        out.update({"c_lambda_own": cl, "c_lambda_project": proj["c_lambda"], "c_simex_own": c_simex,
                    "c_simex_project": proj["c_simex_plan"], "c_simex_absdiff": abs(c_simex - proj["c_simex_plan"]),
                    "simex_x_clipped_own": nclip, "simex_seed": 1110, "Bs": Bs})
        # quadratic extrapolation applied to the PROJECT's c_lambda (deterministic given those values)
        Lp = np.array([float(k) for k in proj["c_lambda"]]); Cp = np.array(list(proj["c_lambda"].values()))
        cp = float(np.polyval(np.polyfit(Lp, Cp, 2), -1.0))
        out["c_simex_from_project_c_lambda"] = cp
        out["c_simex_from_project_c_lambda_absdiff"] = abs(cp - proj["c_simex_plan"])
        out["final_from_project_parts"] = proj["c_simex_plan"] - proj["pb_bias"]
        out["final_from_project_parts_absdiff"] = abs(proj["c_simex_plan"] - proj["pb_bias"] - proj["estimate"])
        print(out, flush=True)
    wjson(out, OUT / f"V11_{EP}_curve.json")


if __name__ == "__main__":
    main()
