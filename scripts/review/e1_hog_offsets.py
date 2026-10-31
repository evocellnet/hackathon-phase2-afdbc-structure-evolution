#!/usr/bin/env python3
"""Review: E1 naive curvature when each HOG gets its own intercept eta_h (y = eta_h + a_h x^c), i.e. allowing a
HOG-specific per-protein 3Di excess (tip noise beyond the identical-sequence floor). Naive fit, no SIMEX.
Also: windows, the eta_h distribution, and a recovery check on noise-free synthetic y = eta + a_h x^0.72.
  CONFIRM_RUN=1 python scripts/review/e1_hog_offsets.py > logs/review/e1_hog_offsets.json"""
import sys, json, tempfile
from pathlib import Path
import numpy as np, polars as pl
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from strucevo.stats.driver import EndpointDriver
from strucevo.stats.endpoints_curve import prepare_curve_data
from scipy.optimize import minimize_scalar

def Q_free(c, g, x, y, w, H):
    u = x ** c
    S = lambda v: np.bincount(g, v, H)
    sw, su, sy, suu, suy = S(w), S(w * u), S(w * y), S(w * u * u), S(w * u * y)
    den = sw * suu - su ** 2
    ok = den > 1e-12 * np.maximum(sw * suu, 1e-300)
    a = np.where(ok, (sw * suy - su * sy) / np.where(ok, den, 1), 0)
    e = np.where(ok, (sy - a * su) / np.maximum(sw, 1e-300), sy / np.maximum(sw, 1e-300))
    r = y - e[g] - a[g] * u
    return float(np.sum(w * r * r)), e, a

def Q_fixed(c, g, x, y, w, H, eta):
    u = x ** c; S = lambda v: np.bincount(g, v, H)
    a = S(w * (y - eta) * u) / np.maximum(S(w * u * u), 1e-300)
    r = y - eta - a[g] * u
    return float(np.sum(w * r * r)), a

def fit(fn):
    grid = np.arange(0.2, 2.01, 0.02)
    q = [fn(c) for c in grid]; i = int(np.argmin(q))
    r = minimize_scalar(fn, bounds=(grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]), method="bounded")
    return float(r.x), float(r.fun), bool(i in (0, len(grid) - 1))

def main():
    d = EndpointDriver(split="confirm", run_alt_alignments=False, outdir=tempfile.mkdtemp(), log=lambda *a: None)
    d.load()
    tab = d.tabs["tdiaa"]; eta = d.tun.eta("tdiaa")["tdi_ml"]
    hc = pl.read_parquet(ROOT / "data/interim/hogs_endpoint_cols_confirm.parquet")
    excl = set(hc.filter(pl.col("excl_tdi_tdiaa"))["hog"].to_list())
    out = {}
    for lo, hi in [(0.05, 2.0), (0.05, 1.0), (0.3, 2.0)]:
        data, info, idx, sp = prepare_curve_data(tab, "E1", excluded_hogs=excl, x_range=(lo, hi))
        m = data.w0 > 0
        g0 = data.g[m]; x = data.x[m]; y = data.y[m]
        ug, g = np.unique(g0, return_inverse=True); H = len(ug)
        w = 1.0 / (data.sy[m] ** 2 + 0.0017371)
        c_fix, q_fix, b1 = fit(lambda c: Q_fixed(c, g, x, y, w, H, eta)[0])
        c_free, q_free, b2 = fit(lambda c: Q_free(c, g, x, y, w, H)[0])
        _, e_h, a_h = Q_free(c_free, g, x, y, w, H)
        _, e1, _ = Q_free(1.0, g, x, y, w, H)
        rec = {"n": int(m.sum()), "H": H, "c_eta_fixed": c_fix, "Q_eta_fixed": q_fix, "c_eta_h_free": c_free,
               "Q_eta_h_free": q_free, "boundary": [b1, b2],
               "Q_eta_h_free_at_c1": Q_free(1.0, g, x, y, w, H)[0],
               "eta_h_at_c1_q10_50_90": np.percentile(e1, [10, 50, 90]).tolist(),
               "eta_h_at_chat_q10_50_90": np.percentile(e_h, [10, 50, 90]).tolist()}
        if (lo, hi) == (0.05, 2.0):
            _, a72 = Q_fixed(0.72, g, x, y, w, H, eta)
            ys = eta + a72[g] * x ** 0.72
            rec["recovery_synthetic_c0.72_eta_h_free"] = fit(lambda c: Q_free(c, g, x, ys, w, H)[0])[0]
            _, a1 = Q_fixed(1.0, g, x, y, w, H, eta)
            for dl in (0.014, 0.03):
                ys2 = eta + 2 * dl + a1[g] * x
                rec[f"synthetic_c1_offset{2*dl:.3f}_eta_fixed_c"] = fit(lambda c: Q_fixed(c, g, x, ys2, w, H, eta)[0])[0]
        out[f"x_{lo}_{hi}"] = rec
        print(json.dumps({f"x_{lo}_{hi}": rec}, default=float), flush=True)

if __name__ == "__main__":
    main()
