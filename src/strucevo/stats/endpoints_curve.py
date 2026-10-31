"""E1 / E2 (ANALYSIS_PLAN §4.1, §4.2): curvature of 3Di-ML (E1) or floor-corrected LDDT distance (E2) vs AA-ML
with the pairwise engine, SIMEX, pigeonhole bootstrap, decision rules, p-value for Holm and d* (E2)."""
from __future__ import annotations

import time

import numpy as np
import torch
from scipy.stats import norm

from .common import (Z98, guard_split, endpoint_record, normal_ci, pct_ci, stage_rng, md5_arrays)
from .pairwise import (CurveData, CurveEngine, two_pass, simex, fit_batch, fit_eta_free, reference_fit, simex_pb_correction,
                       weighted_median, GRID, LAMBDAS)
from .bootstrap import multinomial_weights

CURVE_COLS = {"E1": ("tdi_ml", "tdi_ml_se", "tdi_sat"), "E2": ("lddt_D", "lddt_D_se", "lddt_sat")}


def _get(tab, k):
    if hasattr(tab, "get_column"):
        return tab.get_column(k).to_numpy()
    return np.asarray(tab[k])


def prepare_curve_data(tab, endpoint="E1", excluded_hogs=(), x_range=(0.05, 2.0), min_pairs=10, min_ratio=4.0,
                       ycols=None, require_1to1=True, extra_mask=None, impute_cap=None):
    """Apply the E1/E2 filters. Returns (CurveData, info dict, index arrays)."""
    yc, ysc, ysat = ycols or CURVE_COLS[endpoint]
    guard_split(_get(tab, "split"))
    hog = _get(tab, "hog").astype(str)
    n0 = len(hog)
    keep = np.ones(n0, bool)
    excl = {}
    if require_1to1 and ("one_to_one" in (tab.columns if hasattr(tab, "columns") else tab)):
        m = _get(tab, "one_to_one").astype(bool)
        excl["not_1to1"] = int((keep & ~m).sum()); keep &= m
    xs = _get(tab, "aa_ml_seq").astype(float)
    m = np.isfinite(xs) & (xs >= x_range[0]) & (xs <= x_range[1])
    excl["aa_ml_seq_out_of_range"] = int((keep & ~m).sum()); keep &= m
    m = ~np.isin(hog, np.asarray(list(excluded_hogs), dtype=str)) if len(excluded_hogs) else np.ones(n0, bool)
    excl["hog_excluded_metric"] = int((keep & ~m).sum()); keep &= m
    if extra_mask is not None:
        excl["sensitivity_filter"] = int((keep & ~extra_mask).sum()); keep &= extra_mask
    x = _get(tab, "aa_ml").astype(float)
    y = _get(tab, yc).astype(float)
    sx = _get(tab, "aa_ml_se").astype(float)
    sy = _get(tab, ysc).astype(float)
    status = _get(tab, "status").astype(str)
    sat = _get(tab, ysat).astype(bool) | _get(tab, "aa_sat").astype(bool)
    elig = (status == "ok") & ~sat & np.isfinite(x) & np.isfinite(y) & np.isfinite(sx) & np.isfinite(sy) & (x > 0)
    if impute_cap is not None:       # sensitivity: failed / saturated pairs imputed at the cap
        bad = keep & ~elig
        y = np.where(bad, impute_cap, y)
        sy = np.where(bad & ~np.isfinite(sy), np.nanmedian(sy[elig]), sy)
        x = np.where(bad & ~(np.isfinite(x) & (x > 0)), np.nan, x)
        elig = elig | (bad & np.isfinite(x) & (x > 0) & np.isfinite(sx))
    excl["pair_weight0_fail_or_sat"] = int((keep & ~elig).sum())
    # HOG gate: >= min_pairs eligible pairs and max/min x >= min_ratio
    sel = keep & elig
    hk = hog[sel]
    uh, inv = np.unique(hk, return_inverse=True)
    cnt = np.bincount(inv, minlength=len(uh))
    xmax = np.full(len(uh), -np.inf); xmin = np.full(len(uh), np.inf)
    np.maximum.at(xmax, inv, x[sel]); np.minimum.at(xmin, inv, x[sel])
    good_h = uh[(cnt >= min_pairs) & (xmax / xmin >= min_ratio)]
    m = np.isin(hog, good_h)
    excl["hog_lt_min_pairs_or_range"] = int((keep & ~m).sum())
    keep &= m
    idx = np.flatnonzero(keep)
    hh = hog[idx]
    uh2, g = np.unique(hh, return_inverse=True)
    w0 = elig[idx].astype(float)
    xx = np.where(w0 > 0, x[idx], 1.0)
    yy = np.where(w0 > 0, y[idx], 0.0)
    data = CurveData(g=g, x=xx, y=yy, sx=np.where(w0 > 0, sx[idx], 1.0), sy=np.where(w0 > 0, sy[idx], 1.0),
                     w0=w0, eta=0.0, H=len(uh2), hog_ids=uh2)
    info = {"n_rows_in": n0, "n_pairs_fit": int(w0.sum()), "n_pairs_weight0": int((w0 == 0).sum()),
            "n_hogs": len(uh2), "n_excluded": excl}
    sp = (_get(tab, "species_a").astype(np.int64)[idx], _get(tab, "species_b").astype(np.int64)[idx])
    return data, info, idx, sp


def tree_x(bf, data: CurveData, idx, device=None, metric="aa_ml"):
    """DEVIATION D29: replace x = aa_ml by the tree-smoothed AA distance x_hat = X_h b_hat_AA of the HOG's branch
    fit (pass-2 WLS, calibrated SEs + s2_hm; same fit as E3-E5), and return its error structure for SIMEX.
    idx: row indices (into the pair table the BranchFit was built from) of the curve pairs (prepare_curve_data).
    Pairs of HOGs without an AA fit (underdetermined) or with x_hat <= 0 get weight 0 (counted).
    Modifies `data` in place (x, sx, w0); returns (TreeX, info)."""
    from .pairwise import TreeX, _dev
    dev = _dev(device)
    k = bf.metrics.index(metric)
    by = {hf.hog: hf for hf in bf.hogs}
    n = len(data.x)
    pmax = max([hf.X.shape[1] for hf in bf.hogs] + [1])
    Lp = np.zeros((data.H, pmax, pmax)); Xr = np.zeros((n, pmax))
    xh = np.full(n, np.nan); sd = np.full(n, np.nan); rs = np.zeros(n)
    rowpos = {}
    for hf in bf.hogs:
        for j, r in enumerate(hf.pair_index):
            rowpos[int(r)] = j
    n_nofit = 0
    for h, hid in enumerate(data.hog_ids):
        hf = by.get(str(hid))
        sel = np.flatnonzero(data.g == h)
        if hf is None or hf.fit["under"][k] or not np.all(np.isfinite(hf.fit["b"][k])):
            n_nofit += len(sel); continue
        A = hf.fit["A"][k]; wv = hf.fit["w"][k]; p = A.shape[0]
        okc = wv > 0
        V = (A[:, okc] / wv[okc]) @ A[:, okc].T
        ev_, Q = np.linalg.eigh(0.5 * (V + V.T))
        Lp[h, :p, :p] = Q * np.sqrt(np.clip(ev_, 0, None))
        pos = np.array([rowpos.get(int(r), -1) for r in idx[sel]])
        good = pos >= 0
        Xc = hf.X[pos[good]]
        ss = sel[good]
        Xr[ss, :p] = Xc
        xh[ss] = Xc @ hf.fit["b"][k]
        sd[ss] = np.sqrt(np.einsum("ij,jk,ik->i", Xc, V, Xc))
        own = np.einsum("ij,ji->i", Xc, A[:, pos[good]])                      # hat-matrix own-pair element
        sig_d = np.sqrt(np.where(wv[pos[good]] > 0, 1 / np.where(wv[pos[good]] > 0, wv[pos[good]], 1), 0.0))
        rs[ss] = np.clip(own * sig_d / np.where(sd[ss] > 0, sd[ss], 1), -1, 1)
        n_nofit += int((~good).sum())
    bad = ~(np.isfinite(xh) & (xh > 0) & np.isfinite(sd) & (sd > 0))
    n_nonpos = int((bad & (data.w0 > 0) & np.isfinite(xh)).sum())
    data.w0 = np.where(bad, 0.0, data.w0)
    data.x = np.where(bad, 1.0, xh)
    data.sx = np.where(bad, 1.0, sd)
    t = lambda a: torch.as_tensor(a, dtype=torch.float64, device=dev)
    tx = TreeX(L=t(Lp), Xrow=t(Xr), sd=t(np.where(bad, 1.0, sd)), rho_scale=t(np.where(bad, 0.0, rs)))
    return tx, {"x_tree_n_nofit": int(n_nofit), "x_tree_n_nonpositive": n_nonpos,
                "x_tree_median_sd_ratio": float(np.nanmedian(sd[~bad] / np.maximum(xh[~bad], 1e-9))) if (~bad).any() else None}


def curve_decision(c, se, v3_false_concave=None):
    ci = normal_ci(c, se)
    lo99, hi99 = c - Z98 * se, c + Z98 * se
    if v3_false_concave is not None and v3_false_concave > 0.10:
        cls = "inconclusive"
    elif ci[0] >= 0.9 and ci[1] <= 1.1:
        cls = "linear"
    elif hi99 < 1:
        cls = "concave"
    elif lo99 > 1:
        cls = "inconclusive_convex"
    else:
        cls = "inconclusive"
    p_lin = max(1 - norm.cdf((c - 0.9) / se), 1 - norm.cdf((1.1 - c) / se))
    p_conc = norm.cdf((c - 1) / se)
    p = p_lin if cls == "linear" else p_conc if cls == "concave" else min(p_lin, p_conc)
    return cls, float(p), {"p_lin": float(p_lin), "p_conc": float(p_conc), "lower99": lo99, "upper99": hi99}


def _species_codes(sp_a, sp_b):
    u, inv = np.unique(np.r_[sp_a, sp_b], return_inverse=True)
    return inv[:len(sp_a)], inv[len(sp_a):], len(u)


def fit_curve_full(data: CurveData, eta, kx=1.0, ky=1.0, rho=0.0, B=200, B_s=50, B_s_boot=10, B_naive=1000,
                   seed_stage="bootstrap_E1", sp=None, device=None, window=0.3, simex_on=True, verbose=False,
                   seed_extra=(), xtree=None, pb_K=0, pb_Bs=20, x_design=None):
    """Point estimate (two-pass + SIMEX [+ D29 parametric-bootstrap bias correction]) and pigeonhole bootstrap.
    pb_K > 0: c = c_SIMEX - bias_PB (simex_pb_correction, design x_design or x); bootstrap replicates are shifted
    by the same (fixed) correction and SE = sqrt(SE_boot^2 + MC_SE^2). Returns dict."""
    data.eta = float(eta)
    eng = CurveEngine(data, device=device)
    t0 = time.time()
    f = two_pass(eng)
    out = {"c_naive": f["c"], "s2": f["s2"], "boundary_naive": f["boundary"], "a_h_naive": f["a_h"]}
    rng = stage_rng(seed_stage, *seed_extra, 0)
    if simex_on:
        sm = simex(eng, f["w"], f["c"], kx=kx, ky=ky, rho=rho, B_s=B_s, rng=rng, window=window, xtree=xtree)
        out.update({"c_simex": sm["c_simex"], "c_lambda": sm["c_lambda"], "n_clip": sm["n_clip"],
                    "simex_boundary": sm["boundary_any"]})
        c_hat = sm["c_simex"]
        pb_bias, pb_mc = 0.0, 0.0
        if pb_K:
            pb = simex_pb_correction(eng, f, sm["c_simex"], kx=kx, ky=ky, rho=rho, K=pb_K, B_s=pb_Bs,
                                     rng=stage_rng(seed_stage, *seed_extra, 2), window=window, x_design=x_design)
            pb_bias, pb_mc = pb["bias"], pb["mc_se"]
            out.update({"c_simex_uncorrected": sm["c_simex"], "pb_bias": pb_bias, "pb_mc_se": pb_mc, "pb_K": pb_K,
                        "pb_Bs": pb_Bs, "pb_design": "tree" if x_design is not None else "x"})
            c_hat = pb["c_corrected"]
    else:
        c_hat = f["c"]
    out["c"] = c_hat
    out["a_h"] = eng.a_h(min(max(c_hat, 1e-3), 10), w=f["w"]).cpu().numpy()
    # bootstrap
    g = data.g
    if sp is not None:
        ia, ib, nsp = _species_codes(*sp)
    else:
        ia = ib = np.zeros(len(g), int); nsp = 1
    reps, reps_naive = [], []
    brng = stage_rng(seed_stage, *seed_extra, 1)
    w2 = f["w"]
    for b in range(max(B, B_naive)):
        W = multinomial_weights(brng, data.H)
        M = multinomial_weights(brng, nsp) if sp is not None else np.ones(1)
        pw = W[g] * M[ia] * M[ib]
        wb = w2 * torch.as_tensor(pw, dtype=torch.float64, device=eng.dev)
        if b < B_naive:
            cn, _ = fit_batch(eng, w=wb, window=(f["c"] - window, f["c"] + window))
            reps_naive.append(cn[0])
        if b < B and simex_on:
            sm = simex(eng, wb, cn[0] if b < B_naive else fit_batch(eng, w=wb)[0][0], kx=kx, ky=ky, rho=rho,
                       B_s=B_s_boot, rng=brng, window=window, xtree=xtree)
            reps.append(sm["c_simex"])
        if verbose and b % 50 == 0:
            print(f"   boot {b}  {time.time() - t0:.1f}s", flush=True)
    reps = np.asarray(reps if simex_on else reps_naive[:B])
    if simex_on and pb_K:
        reps = reps - pb_bias                      # fixed correction (held, like the pass-2 weights, D5)
    out["reps"] = reps
    out["reps_naive"] = np.asarray(reps_naive)
    out["se_boot"] = float(np.sqrt(np.std(reps, ddof=1) ** 2 + (pb_mc ** 2 if simex_on and pb_K else 0.0)))
    out["se_boot_naive"] = float(np.std(reps_naive, ddof=1))
    out["runtime_s"] = time.time() - t0
    out["engine"] = eng
    out["w"] = w2
    return out


def d_star(x, y_raw, w, g, H, eta_raw, x_grid=np.round(np.arange(0.05, 2.0 + 1e-9, 0.01), 2), engine=None):
    """E2 d*: profiled linear fit eta + a_h x on x <= 0.5 and power fit eta + a'_h x^c' on x <= 2 for
    y_raw = 1 - LDDT_global; A, A' = weighted medians of a_h, a'_h (weights = per-HOG Fisher information
    sum w x^2 resp. sum w x^{2c'}; clarification C4); d* = smallest grid x with |power/linear - 1| > 0.10."""
    lin = x <= 0.5
    num = np.bincount(g[lin], w[lin] * (y_raw[lin] - eta_raw) * x[lin], H)
    den = np.bincount(g[lin], w[lin] * x[lin] ** 2, H)
    a_lin = np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)
    A = weighted_median(a_lin, den)
    pw = x <= 2.0
    if engine is None:
        cd = CurveData(g=g[pw], x=x[pw], y=y_raw[pw], sx=np.ones(pw.sum()), sy=np.ones(pw.sum()), w0=np.ones(pw.sum()),
                       eta=eta_raw, H=H)
        engine = CurveEngine(cd)
    wt = torch.as_tensor(w[pw], dtype=torch.float64, device=engine.dev)
    cp = fit_batch(engine, w=wt[None])[0][0]
    xc = x[pw] ** cp
    num2 = np.bincount(g[pw], w[pw] * (y_raw[pw] - eta_raw) * xc, H)
    den2 = np.bincount(g[pw], w[pw] * xc * xc, H)
    a_pow = np.where(den2 > 0, num2 / np.where(den2 > 0, den2, 1), np.nan)
    Ap = weighted_median(a_pow, den2)
    rel = np.abs((eta_raw + Ap * x_grid ** cp) / (eta_raw + A * x_grid) - 1)
    hit = np.flatnonzero(rel > 0.10)
    ok = np.flatnonzero(rel <= 0.10)
    # descriptive only (not the plan's estimand): where the two fits agree within 10%, and the departure above it
    desc = {"rel_at_grid_min": float(rel[0]),
            "agree_within_10pct": [float(x_grid[ok[0]]), float(x_grid[ok[-1]])] if ok.size else None,
            "upper_departure_x": (float(x_grid[ok[-1] + 1]) if ok.size and ok[-1] + 1 < len(x_grid) else None)}
    return {"d_star": float(x_grid[hit[0]]) if hit.size else float("nan"), "A": A, "A_pow": Ap, "c_pow": float(cp),
            "at_grid_min": bool(hit.size and hit[0] == 0), "descriptive": desc}, engine


def d_star_boot(x, y_raw, w, g, H, eta_raw, sp, B=200, seed_stage="bootstrap_E2"):
    pt, eng = d_star(x, y_raw, w, g, H, eta_raw)
    ia, ib, nsp = _species_codes(*sp)
    rng = stage_rng(seed_stage, 7)
    reps = []
    for _ in range(B):
        W = multinomial_weights(rng, H); M = multinomial_weights(rng, nsp)
        reps.append(d_star(x, y_raw, w * W[g] * M[ia] * M[ib], g, H, eta_raw, engine=eng)[0]["d_star"])
    reps = np.asarray(reps)
    se = float(np.nanstd(reps, ddof=1))
    pt["se_boot"] = se
    pt["ci98"] = normal_ci(pt["d_star"], se)
    pt["ci_pct98"] = pct_ci(reps)
    pt["n_rep_nan"] = int(np.sum(~np.isfinite(reps)))
    return pt


def run_curve_endpoint(tab, endpoint="E1", eta=0.0, kx=1.0, ky=1.0, rho=0.0, excluded_hogs=(), B=200, B_s=50,
                       B_s_boot=10, B_naive=1000, v3_false_concave=None, sensitivities=None, B_sens=50,
                       eta_raw=None, device=None, alignment_results=None, verbose=False, input_md5=None,
                       bf=None, rho_pair=None, pb_K=25, pb_Bs=20, pb_K_sens=10):
    """Full E1/E2 endpoint. `sensitivities` = {name: dict(tab=..., extra_mask=..., eta_free=True, impute_cap=...,
    x_range=..., ycols=..., bf=...)}; each is re-estimated (SIMEX + B_sens bootstrap) and reported with its decision.
    `alignment_results` = {aln: (estimate, se)} for seqonly / tdiaa (alignment robustness).
    DEVIATION D29: the SIMEX estimate is bias-corrected by a parametric bootstrap (pb_K pseudo worlds,
    simex_pb_correction); the plan's uncorrected SIMEX value is reported as `c_simex_plan`. bf: BranchFit built from
    THIS `tab` (row order); used ONLY for the pseudo-world design D* = tree-smoothed AA distance (tree_x); without
    it D* = observed x. The regressor is always x = aa_ml (plan). rho_pair: per-row rho(bin) array over `tab`
    (per-pair Sigma_ij); None -> scalar `rho`. pb_K = 0 reproduces the plan's estimator."""
    stage = "bootstrap_E1" if endpoint == "E1" else "bootstrap_E2"
    import copy

    def _prep(t, bfx, rp, **kw):
        d, inf, ix, spx = prepare_curve_data(t, endpoint, **kw)
        xd = None
        if bfx is not None:
            dc = copy.deepcopy(d)
            _, xi = tree_x(bfx, dc, ix, device=device)
            xd = np.where(dc.w0 > 0, dc.x, d.x)
            inf["pb_design_tree"] = xi
        r_ = rho if rp is None else np.asarray(rp, float)[ix]
        return d, inf, ix, spx, xd, r_

    data, info, idx, sp, xdes, rho_used = _prep(tab, bf, rho_pair, excluded_hogs=excluded_hogs)
    res = fit_curve_full(data, eta, kx, ky, rho_used, B=B, B_s=B_s, B_s_boot=B_s_boot, B_naive=B_naive,
                         seed_stage=stage, sp=sp, device=device, verbose=verbose, pb_K=pb_K, pb_Bs=pb_Bs,
                         x_design=xdes)
    xtree = None
    c, se = res["c"], res["se_boot"]
    cls, p, extra = curve_decision(c, se, v3_false_concave)
    flags = []
    if res["boundary_naive"]:
        flags.append("grid_min_at_boundary")
    if res.get("n_clip", 0):
        flags.append(f"simex_x_clipped:{res['n_clip']}")
    ci_n = normal_ci(c, se)
    ci_p = pct_ci(res["reps"])
    cls_p, _, _ = curve_decision(c, (ci_p[1] - ci_p[0]) / (2 * Z98) if np.all(np.isfinite(ci_p)) else se, v3_false_concave)
    if cls_p != cls:
        flags.append("ci_disagree")
    # independent recomputation (CPU numpy) of the naive estimate
    w = res["w"].cpu().numpy()
    c_ref = reference_fit(data.g, data.x, data.y, w, data.eta)
    check_recompute = abs(c_ref - res["c_naive"])
    if check_recompute > 1e-6:
        flags.append("recompute_mismatch")
    sens = {}
    sens_items = list((sensitivities or {}).items())
    for name, opt in sens_items:
        opt = dict(opt)
        note = opt.pop("note", None)
        same_tab = "tab" not in opt
        t2 = opt.pop("tab", tab)
        bf2 = opt.pop("bf", bf if same_tab else None)
        if opt.pop("eta_free", False):
            f2 = fit_eta_free(res["engine"], res["w"])
            sens[name] = {"estimate": f2["c"], "eta": f2["eta"], "ci98": None, "decision": None,
                          "note": "naive c with eta profiled; no bootstrap"}
            continue
        rp2 = rho_pair if same_tab else None
        r_sc = rho if rho_pair is None else float(np.mean(np.asarray(rho_pair, float)[idx]))
        d2, i2, _, sp2, xd2, rho2 = _prep(t2, bf2, rp2, excluded_hogs=opt.pop("excluded_hogs", excluded_hogs), **opt)
        if rp2 is None:
            rho2 = r_sc
        r2 = fit_curve_full(d2, eta, kx, ky, rho2, B=B_sens, B_s=B_s, B_s_boot=B_s_boot, B_naive=B_sens,
                            seed_stage=stage, sp=sp2, device=device, pb_K=pb_K_sens if pb_K else 0, pb_Bs=pb_Bs,
                            x_design=xd2)
        c2, s2 = r2["c"], r2["se_boot"]
        cl2, _, _ = curve_decision(c2, s2, v3_false_concave)
        sens[name] = {"estimate": c2, "ci98": normal_ci(c2, s2), "decision": cl2, "n_pairs": i2["n_pairs_fit"],
                      "c_simex_plan": r2.get("c_simex_uncorrected"), "pb_bias": r2.get("pb_bias")}
        if note:
            sens[name]["note"] = note
    arob = None
    if alignment_results:
        classes = [cls] + [curve_decision(e, s, v3_false_concave)[0] for e, s in alignment_results.values()]
        signs = [np.sign(c - 1)] + [np.sign(e - 1) for e, s in alignment_results.values()]
        arob = bool(len(set(classes)) == 1 and len(set(signs)) == 1)
        for a_, (e, s) in alignment_results.items():
            sens[f"aln_{a_}"] = {"estimate": e, "ci98": normal_ci(e, s), "decision": curve_decision(e, s, v3_false_concave)[0]}
    rec = endpoint_record(
        endpoint=endpoint,
        estimand="global curvature c in y = eta + a_h x^c (SIMEX-corrected)" + (" ; y = tdi_ml" if endpoint == "E1" else " ; y = lddt_D")
                 + (" ; SIMEX bias-corrected by parametric bootstrap (D29)" if pb_K else ""),
        estimate=c, se_boot=se, ci98=ci_n, ci_pct98=ci_p,
        one_sided_bound={"lower99": extra["lower99"], "upper99": extra["upper99"]},
        p_value=p, p_holm=None, decision=cls, n_hogs=info["n_hogs"], n_units=info["n_pairs_fit"],
        n_excluded=info["n_excluded"], sensitivity=sens, alignment_robust=arob, flags=flags,
        input_md5=input_md5 or md5_arrays(data.g, data.x, data.y),
        c_naive=res["c_naive"], c_naive_se_boot_B1000=res["se_boot_naive"], c_lambda=res.get("c_lambda"),
        s2=res["s2"], p_lin=extra["p_lin"], p_conc=extra["p_conc"], v3_false_concave=v3_false_concave,
        checks={"recompute_abs_diff": check_recompute, "grid_boundary": res["boundary_naive"],
                "n_fit_equals_filter_log": info["n_pairs_fit"] == int(data.w0.sum())},
        eta=eta, kappa_x=kx, kappa_y=ky, rho=rho if rho_pair is None else "per_pair_bin",
        rho_mean=float(np.mean(rho_used)), c_simex_plan=res.get("c_simex_uncorrected"), pb_bias=res.get("pb_bias"),
        pb_mc_se=res.get("pb_mc_se"), pb_K=pb_K, pb_design=res.get("pb_design"), pb_design_tree=info.get("pb_design_tree"),
        B=B, B_s=B_s, B_s_boot=B_s_boot, runtime_s=res["runtime_s"])
    if endpoint == "E2" and eta_raw is not None:
        lg = _get(tab, "lddt_global").astype(float)[idx]
        rec["d_star"] = d_star_boot(data.x, 1 - lg, (data.w0 > 0).astype(float) / data.sy ** 2, data.g, data.H,
                                    eta_raw, sp, B=B, seed_stage=stage)
    curve = {"c": c, "eta": eta, "a_h": dict(zip(data.hog_ids.tolist(), res["a_h"].tolist()))}
    return rec, curve, res


def standard_sensitivities(tab, hidden_para_rows=None, single_copy_hogs=None, variants=None):
    """Plan §4.1 sensitivity set for E1/E2. `variants` = {name: alternative pair table} for the ones that need
    recomputed inputs (all aligned columns, alpha x1.5 / /1.5); pass `alignment_results` separately."""
    n = len(_get(tab, "hog"))
    s = {"eta_free": {"eta_free": True},
         "d_aa_le_1": {"x_range": (0.05, 1.0)},
         "failures_imputed_at_cap": {"impute_cap": 20.0}}
    if hidden_para_rows is not None:
        m = np.ones(n, bool); m[np.asarray(hidden_para_rows, int)] = False
        s["hidden_paralogy_removed"] = {"extra_mask": m}
    if single_copy_hogs is not None:
        s["single_copy_hogs"] = {"extra_mask": np.isin(_get(tab, "hog").astype(str), np.asarray(list(single_copy_hogs), dtype=str))}
    for k, v in (variants or {}).items():
        s[k] = {"tab": v}
    return s
