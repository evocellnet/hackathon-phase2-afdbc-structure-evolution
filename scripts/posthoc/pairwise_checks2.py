#!/usr/bin/env python3
"""Post-hoc verification, part 2: errors-in-variables matter for every naive estimator, so
 (a) real-data curvature fits on the tree-smoothed AA distance x_hat = X_h b_hat_AA (production D29 x), three baselines;
 (b) synthetic recovery on the real design with REALISTIC, mostly tree-additive x errors:
     x_true = x_hat (tree-additive); x_obs = x_true + X_h e (e ~ N(0, V_b) from the production AA branch fit)
     [+ iid remainder with 17% of the pair-level Fisher variance, mimicking the measured 13-21% non-additive share];
     y = eta + u_i + u_j + a_h x_true^c + N(0, sy^2 + s2); u_i = 0 or Gamma(1, 0.015) (mean 0.015 per protein).
  CONFIRM_RUN=1 python scripts/posthoc/pairwise_checks2.py -> results/posthoc/pairwise_checks2.json
"""
import sys, json, time, tempfile
from pathlib import Path
import numpy as np, polars as pl
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(Path(__file__).resolve().parent))
from strucevo.stats.driver import EndpointDriver
from common import get_window, Fitter

OUT = ROOT / "results" / "posthoc"
t0 = time.time()
def log(*a): print(f"[{time.time()-t0:7.0f}s]", *a, flush=True)
res = {}
def save(): (OUT / "pairwise_checks2.json").write_text(json.dumps(res, indent=1, default=float))

d = EndpointDriver(split="confirm", run_alt_alignments=False, outdir=tempfile.mkdtemp(), log=lambda *a: None)
d.load(); d.tau = d.tun.tau_misfit()
tab = d.tabs["tdiaa"]; ETA = d.tun.eta("tdiaa")["tdi_ml"]
hc = pl.read_parquet(ROOT / "data/interim/hogs_endpoint_cols_confirm.parquet")
excl = set(hc.filter(pl.col("excl_tdi_tdiaa"))["hog"].to_list())
bf = d.fit(tab, "tdiaa")
log("branch fit done")
rowpos = {}
for hf in bf.hogs:
    for j, r in enumerate(hf.pair_index):
        rowpos[int(r)] = (hf, j)

def smooth(D):
    """x_hat and per-HOG (X rows, covariance factor) for the curve pairs."""
    n = len(D["x"]); xh = np.full(n, np.nan); info = {}
    hogs_rows = {}
    for i, r in enumerate(D["rows"]):
        t = rowpos.get(int(r))
        if t is None: continue
        hf, j = t
        if hf.fit["under"][0] or not np.all(np.isfinite(hf.fit["b"][0])): continue
        hogs_rows.setdefault(id(hf), (hf, []))[1].append((i, j))
    for hf, lst in hogs_rows.values():
        ii = np.array([a for a, _ in lst]); jj = np.array([b for _, b in lst])
        Xc = hf.X[jj]; b = hf.fit["b"][0]
        A = hf.fit["A"][0]; wv = hf.fit["w"][0]; okc = wv > 0
        V = (A[:, okc] / wv[okc]) @ A[:, okc].T
        ev, Q = np.linalg.eigh(0.5 * (V + V.T)); L = Q * np.sqrt(np.clip(ev, 0, None))
        xh[ii] = Xc @ b
        info[hf.hog] = (ii, Xc, L)
    return xh, info

res["fits_xhat"] = {}
for lo, hi in [(0.05, 2.0), (0.3, 2.0)]:
    D = get_window(tab, excl, lo, hi)
    xh, info = smooth(D)
    keep = np.isfinite(xh) & (xh > 1e-3)
    # restrict D to pairs with a usable x_hat
    sel = np.flatnonzero(keep)
    D2 = {k: (v[sel] if isinstance(v, np.ndarray) and len(v) == len(keep) else v) for k, v in D.items()}
    D2["xx"] = xh[sel]
    f0 = Fitter(D2, 1 / D2["sy"] ** 2, "fixed"); c0 = f0.fit_c(D2["y"], ETA)[0]
    r0 = f0.Q(c0, D2["y"], ETA)[2]; s2 = max(0.0, float(np.mean(r0 ** 2) - np.mean(D2["sy"] ** 2)))
    w = 1 / (D2["sy"] ** 2 + s2)
    rec = {"n": len(sel), "n_dropped_no_xhat": int((~keep).sum()), "s2": s2}
    Fs = {}
    for mode in ("fixed", "hog", "prot"):
        F = Fitter(D2, w, mode); Fs[mode] = F
        c, q, b = F.fit_c(D2["y"], ETA)
        rec[mode] = {"c": c, "Q": q, "boundary": b, "Q_c1_minus_Qmin": F.Q(1.0, D2["y"], ETA)[0] - q}
        if mode != "fixed":
            off = F.offsets(c, D2["y"], ETA); off1 = F.offsets(1.0, D2["y"], ETA)
            rec[mode]["offset_median_at_chat"] = float(np.median(off)); rec[mode]["offset_median_at_c1"] = float(np.median(off1))
        log("xhat", lo, hi, mode, json.dumps(rec[mode], default=float))
    res["fits_xhat"][f"x_{lo}_{hi}"] = rec; save()
    if (lo, hi) != (0.05, 2.0):
        continue
    # ------------------------------------------------ synthetic on this design
    # map info (indices into D) to indices into D2
    newpos = np.full(len(keep), -1); newpos[sel] = np.arange(len(sel))
    blocks = []
    for h, (ii, Xc, L) in info.items():
        m = keep[ii]
        if m.any():
            blocks.append((newpos[ii[m]], Xc[m], L))
    x_true = D2["xx"].copy(); sx = D2["sx"]; n = len(x_true); P = D2["P"]
    rng = np.random.default_rng(31102026)
    syn = []
    for c_true in (0.72, 1.0):
        _, a_true, _ = Fs["fixed"].Q(c_true, D2["y"], ETA)
        for off_kind in ("none", "prot_gamma_mean0.015"):
            for noise in ("tree_additive", "tree_additive+iid17"):
                for rep in range(2):
                    uu = np.zeros(P) if off_kind == "none" else rng.gamma(1.0, 0.015, P)
                    ys = ETA + uu[D2["pa"]] + uu[D2["pb"]] + a_true[D2["g"]] * x_true ** c_true \
                        + rng.normal(0, 1, n) * np.sqrt(D2["sy"] ** 2 + s2)
                    ex = np.zeros(n)
                    for ii, Xc, L in blocks:
                        ex[ii] = Xc @ (L @ rng.normal(0, 1, L.shape[1]))
                    if noise.endswith("iid17"):
                        ex += rng.normal(0, 1, n) * sx * np.sqrt(0.17)
                    D2["xx"] = np.maximum(x_true + ex, 1e-4)
                    r = {"c_true": c_true, "offsets": off_kind, "noise": noise, "rep": rep,
                         "x_err_sd_median": float(np.median(np.abs(ex))) * 1.4826,
                         "c_fixed": Fs["fixed"].fit_c(ys, ETA)[0], "c_hog": Fs["hog"].fit_c(ys)[0],
                         "c_prot": Fs["prot"].fit_c(ys)[0]}
                    D2["xx"] = x_true
                    syn.append(r); log("syn", json.dumps(r)); res["synthetic_tree_noise"] = syn; save()
log("done")
