"""E3, E4, E5 (ANALYSIS_PLAN §4.3-§4.5) on the branch-engine output (BranchFit)."""
from __future__ import annotations

import time
import warnings
from dataclasses import dataclass

import numpy as np
from scipy.stats import norm

from .common import Z98, endpoint_record, normal_ci, pct_ci, stage_rng, md5_arrays
from .bootstrap import multinomial_weights, replicate_groups
from .trees import PanelTree, draw_ages
from .reml import fast_mixed

BIN_EDGES = np.array([0.0, 30.0, 100.0, 300.0, 1000.0, np.inf])


# ------------------------------------------------------------------------------------------ edge view
@dataclass
class EdgeView:
    """AA and 3Di edge rows aligned (one row per HOG edge), restricted to included HOGs."""
    hog_code: np.ndarray
    hog_ids: np.ndarray            # (H,) ids for codes
    child: np.ndarray              # panel node index of child
    parent: np.ndarray
    is_terminal: np.ndarray
    simple: np.ndarray             # maps to exactly one panel edge (and not a merged root edge)
    b_aa: np.ndarray
    b_t: np.ndarray
    se_aa: np.ndarray
    se_t: np.ndarray
    cov: np.ndarray
    span_ptr: np.ndarray           # CSR of spanned panel edges (child-node ids)
    span_edges: np.ndarray
    hog_var_aa: np.ndarray         # (H,) var of HOG total
    hog_var_t: np.ndarray
    hog_cov: np.ndarray
    n_excluded: dict


def edge_view(bf, panel: PanelTree, metric_t="tdi_ml", tau_misfit=None, hog_subset=None, internal_only=False,
              drop_hogs=()):
    E = bf.edges
    F = bf.fits
    ma = E["metric"] == "aa_ml"
    mt = E["metric"] == metric_t
    # HOG-level exclusion: any reason in aa or the 3Di metric, or misfit
    excl_h = {}
    for m in ("aa_ml", metric_t):
        s = F["metric"] == m
        for h, ex, rs, mr in zip(F["hog"][s], F["excluded"][s], F["reason"][s], F["misfit_ratio"][s]):
            reasons = [r for r in str(rs).split(";") if r]
            if tau_misfit is not None and np.isfinite(mr) and mr > tau_misfit and "misfit" not in reasons:
                reasons.append("misfit")
            if reasons:
                excl_h.setdefault(h, set()).update(reasons)
    for h in drop_hogs:
        excl_h.setdefault(h, set()).add("sensitivity_drop")
    hogs_a = E["hog"][ma]
    keep = ~np.isin(hogs_a, np.array(list(excl_h.keys()), dtype=object)) if excl_h else np.ones(ma.sum(), bool)
    if hog_subset is not None:
        keep &= np.isin(hogs_a, np.asarray(list(hog_subset), dtype=object))
    if internal_only:
        keep &= ~E["is_terminal"][ma]
    n_excl = {}
    for h, rs in excl_h.items():
        for r in rs:
            n_excl[r] = n_excl.get(r, 0) + 1
    ia = np.flatnonzero(ma)[keep]
    it = np.flatnonzero(mt)[keep]
    assert np.all(E["hog"][ia] == E["hog"][it]) and np.all(E["_child_idx"][ia] == E["_child_idx"][it])
    hog_ids, code = np.unique(E["hog"][ia], return_inverse=True)
    child = E["_child_idx"][ia].astype(int)
    parent = E["_parent_idx"][ia].astype(int)
    # spans
    ptr = [0]; spans = []
    for c, p in zip(child, parent):
        while c != p:
            spans.append(c); c = panel.parent[c]
        ptr.append(len(spans))
    # HOG total variances
    fa = (F["metric"] == "aa_ml"); ft = (F["metric"] == metric_t)
    va = dict(zip(F["hog"][fa], F["var_sum"][fa]))
    vt = dict(zip(F["hog"][ft], F["var_sum"][ft]))
    cv = dict(zip(F["hog"][ft], F["cov_sum_aa"][ft]))
    simple = np.array([pe is not None for pe in E["panel_edge"][ia]])
    return EdgeView(hog_code=code, hog_ids=hog_ids, child=child, parent=parent, is_terminal=E["is_terminal"][ia],
                    simple=simple, b_aa=E["b"][ia].astype(float), b_t=E["b"][it].astype(float),
                    se_aa=E["b_se"][ia].astype(float), se_t=E["b_se"][it].astype(float), cov=E["cov_aa"][it].astype(float),
                    span_ptr=np.asarray(ptr), span_edges=np.asarray(spans, dtype=int),
                    hog_var_aa=np.array([va[h] for h in hog_ids]), hog_var_t=np.array([vt[h] for h in hog_ids]),
                    hog_cov=np.array([cv[h] for h in hog_ids]), n_excluded=n_excl)


def _edge_w(ev: EdgeView, W, L):
    lam = L[ev.span_edges]
    s = np.add.reduceat(lam, ev.span_ptr[:-1]) if len(lam) else np.zeros(len(ev.child))
    cnt = np.diff(ev.span_ptr)
    s = np.where(cnt > 0, s, 0.0)
    return W[ev.hog_code] * s / np.maximum(cnt, 1)


# ------------------------------------------------------------------------------------------ E3
def _e3_stats(ev, w, hog_w, cls_code=None, groups=None):
    H = len(ev.hog_ids)
    SA = np.bincount(ev.hog_code, w * ev.b_aa, H)
    ST = np.bincount(ev.hog_code, w * ev.b_t, H)
    k = SA.sum() / ST.sum()
    # per-HOG k_h (unweighted sums within the HOG -- W_h enters as a HOG weight)
    SA0 = SA / np.maximum(hog_w, 1e-300)          # per-HOG Lambda-weighted sums (W_h divided out)
    ST0 = ST / np.maximum(hog_w, 1e-300)
    ok = (SA0 > 0) & (ST0 > 0) & (hog_w > 0)
    lk = np.log(SA0[ok] / ST0[ok])
    v = (ev.hog_var_aa / SA0 ** 2 + ev.hog_var_t / ST0 ** 2 - 2 * ev.hog_cov / (SA0 * ST0))[ok]
    hw = hog_w[ok]
    mu = np.average(lk, weights=hw)
    var = np.average((lk - mu) ** 2, weights=hw) * hw.sum() / max(hw.sum() - 1, 1)
    s2 = var - np.average(v, weights=hw)
    out = {"k": k, "log_k": np.log(k), "sigma2_k": s2, "sigma_k": np.sqrt(max(0.0, s2)), "n_hog_logk": int(ok.sum())}
    if cls_code is not None:
        for name, (a, b) in {"cath_alpha_vs_beta": ("alpha", "beta")}.items():
            ma = cls_code == a; mb = cls_code == b
            ka = SA[ma].sum() / ST[ma].sum() if ma.any() else np.nan
            kb = SA[mb].sum() / ST[mb].sum() if mb.any() else np.nan
            out["k_alpha"], out["k_beta"] = ka, kb
            out["delta_cath"] = np.log(ka) - np.log(kb) if (ka > 0 and kb > 0) else np.nan
    for gname, gcode in (groups or {}).items():
        for lev in np.unique(gcode):
            m = gcode == lev
            out[f"k_{gname}_{lev}"] = SA[m].sum() / ST[m].sum()
    return out


def run_E3(bf, panel, tau_misfit=None, cath_class=None, hog_groups=None, B=1000, alignment_k=None, input_md5=None):
    """E3: k = sum b_AA / sum b_3Di (tip-corrected), sigma_k, CATH contrast, subgroup k's, internal-only k."""
    t0 = time.time()
    ev = edge_view(bf, panel, tau_misfit=tau_misfit)
    H = len(ev.hog_ids)
    cls = np.array([cath_class.get(h, "none") for h in ev.hog_ids]) if cath_class else None
    groups = {g: np.array([d.get(h, "NA") for h in ev.hog_ids]) for g, d in (hog_groups or {}).items()}
    ones_w = np.ones(len(ev.child))
    pt = _e3_stats(ev, ones_w, np.ones(H), cls, groups)
    ev_int = edge_view(bf, panel, tau_misfit=tau_misfit, internal_only=True)
    pt_int = _e3_stats(ev_int, np.ones(len(ev_int.child)), np.ones(len(ev_int.hog_ids)))
    rng = stage_rng("bootstrap_E3", 0)
    nE = panel.n_nodes
    reps = []; reps_int = []
    for b in range(B):
        W = multinomial_weights(rng, H)
        L = np.zeros(nE); pe = panel.panel_edge_ids(); L[pe] = multinomial_weights(rng, len(pe))
        reps.append(_e3_stats(ev, _edge_w(ev, W, L), W, cls, groups))
        Wi = np.zeros(len(ev_int.hog_ids))
        # map W onto the internal-only view's HOG codes
        pos = {h: i for i, h in enumerate(ev.hog_ids)}
        Wi = np.array([W[pos[h]] if h in pos else 0.0 for h in ev_int.hog_ids])
        reps_int.append(_e3_stats(ev_int, _edge_w(ev_int, Wi, L), Wi)["k"])
    R = {key: np.array([r.get(key, np.nan) for r in reps]) for key in pt}
    # CI on log k, back-transformed (k is a ratio)
    se_lk = np.std(np.log(R["k"]), ddof=1)
    ci_k = list(np.exp(normal_ci(np.log(pt["k"]), se_lk)))
    ci_pct = pct_ci(R["k"])
    out_sens = {}
    for key in pt:
        if key.startswith("k_") and key not in ("k_alpha", "k_beta"):
            s = np.std(np.log(R[key][R[key] > 0]), ddof=1) if np.sum(R[key] > 0) > 2 else np.nan
            out_sens[key] = {"estimate": pt[key], "ci98": list(np.exp(normal_ci(np.log(pt[key]), s))) if pt[key] > 0 else None, "decision": None}
    sik = np.std(np.log(np.array(reps_int)), ddof=1)
    out_sens["internal_edges_only"] = {"estimate": pt_int["k"], "ci98": list(np.exp(normal_ci(np.log(pt_int["k"]), sik))), "decision": None}
    cath = None
    if cls is not None and np.isfinite(pt.get("delta_cath", np.nan)):
        sd = np.nanstd(R["delta_cath"], ddof=1)
        z = pt["delta_cath"] / sd
        cath = {"delta": pt["delta_cath"], "se_boot": sd, "ci95": normal_ci(pt["delta_cath"], sd, z=1.959964),
                "p_two_sided_unadjusted": float(2 * (1 - norm.cdf(abs(z)))), "k_alpha": pt["k_alpha"], "k_beta": pt["k_beta"],
                "n_alpha": int((cls == "alpha").sum()), "n_beta": int((cls == "beta").sum())}
    sig_se = np.nanstd(R["sigma_k"], ddof=1)
    arob = None
    if alignment_k:
        arob = bool(all(abs(v / pt["k"] - 1) <= 0.25 for v in alignment_k.values()))
        for a_, v in alignment_k.items():
            out_sens[f"aln_{a_}"] = {"estimate": v, "ci98": None, "decision": None}
    rec = endpoint_record(
        endpoint="E3", estimand="k = sum b_AA / sum b_3Di (tip-corrected, all edges, ratio of sums)",
        estimate=pt["k"], se_boot=float(np.std(R["k"], ddof=1)), ci98=ci_k, ci_pct98=ci_pct,
        one_sided_bound=None, p_value=None, p_holm=None, decision="estimation_only",
        n_hogs=H, n_units=int(len(ev.child)), n_excluded=ev.n_excluded, sensitivity=out_sens,
        alignment_robust=arob, flags=[], input_md5=input_md5 or md5_arrays(ev.b_aa, ev.b_t),
        sigma_k=pt["sigma_k"], sigma2_k=pt["sigma2_k"], sigma_k_ci98=normal_ci(pt["sigma_k"], sig_se),
        ci98_overlaps_3_10=bool(ci_k[1] >= 3 and ci_k[0] <= 10), cath_contrast=cath, B=B, runtime_s=time.time() - t0,
        r_ratio=1 / pt["k"])
    return rec


# ------------------------------------------------------------------------------------------ E4
def build_cells(ev: EdgeView, age: np.ndarray, w=None, min_aa=0.05, drop_bin5=False, tau_perm=None,
                rule="aa", info_T=None, tip_var=0.0, r_bar=None, tip_all=False):
    """Cells (h, g) with the AA-only merge rule (ANALYSIS_PLAN §2.6). age: per panel node.

    rule = 'aa'   : the plan's rule (merge while S_AA < min_aa).
    rule = 'info' : DEVIATION D30 (agents/bias_fix_report.md): additionally merge while the cell's PREDICTED 3Di
                    signal-to-noise r_bar * S_AA / sqrt(V_t) < info_T, where V_t = sum w^2 (b_se_3Di^2 + tip_var *
                    terminal) (3Di SEs + per-protein tip-noise variance; no 3Di VALUES of the cell enter) and r_bar
                    is one global scalar (the D8 precedent). Merge direction as the plan (into the next-older cell;
                    the oldest into the next-younger).
    Returns dict of arrays: hog, S_aa, S_t, n_edges, mean_tau, censored, bin, V_aa, V_t (delta-method variances of
    the cell sums: sum w^2 se^2, plus tip_var on terminal 3Di edges)."""
    w = np.ones(len(ev.child)) if w is None else w
    tau = 0.5 * (age[ev.child] + age[ev.parent])
    if tau_perm is not None:          # V6 'edges' negative control: tau permuted across all HOG edges
        tau = tau[tau_perm]
    g = np.clip(np.searchsorted(BIN_EDGES, tau, side="right") - 1, 0, 4)
    H = len(ev.hog_ids)
    key = ev.hog_code * 5 + g
    SA = np.bincount(key, w * ev.b_aa, H * 5).reshape(H, 5)
    ST = np.bincount(key, w * ev.b_t, H * 5).reshape(H, 5)
    NE = np.bincount(key, w, H * 5).reshape(H, 5)
    TS = np.bincount(key, w * tau, H * 5).reshape(H, 5)
    se_a = np.nan_to_num(ev.se_aa, nan=0.0); se_t = np.nan_to_num(ev.se_t, nan=0.0)
    VA = np.bincount(key, w * w * se_a ** 2, H * 5).reshape(H, 5)
    tip_w = np.ones(len(ev.child)) if tip_all else ev.is_terminal.astype(float)
    VT = np.bincount(key, w * w * (se_t ** 2 + tip_var * tip_w), H * 5).reshape(H, 5)
    present = np.bincount(key, (w > 0).astype(float), H * 5).reshape(H, 5) > 0
    out = {k: [] for k in ("hog", "S_aa", "S_t", "n_edges", "mean_tau", "bin", "V_aa", "V_t")}
    n_aa_nonpos = 0
    if rule == "info" and r_bar is None:
        r_bar = float(np.sum(w * ev.b_t) / np.sum(w * ev.b_aa))
    for h in np.flatnonzero(present.any(1)):
        cells = [[j, SA[h, j], ST[h, j], NE[h, j], TS[h, j], VA[h, j], VT[h, j]] for j in range(5) if present[h, j]]
        while len(cells) > 1:
            if rule == "info":
                bad = [i for i, c in enumerate(cells)
                       if c[1] < min_aa or r_bar * c[1] < info_T * np.sqrt(max(c[6], 0.0))]
            else:
                bad = [i for i, c in enumerate(cells) if c[1] < min_aa]
            if not bad:
                break
            i = bad[0]
            tgt = i + 1 if i + 1 < len(cells) else i - 1
            for q in (1, 2, 3, 4, 5, 6):
                cells[tgt][q] += cells[i][q]
            del cells[i]
        for c in cells:
            if c[1] <= 0:
                n_aa_nonpos += 1
                continue
            out["hog"].append(h); out["S_aa"].append(c[1]); out["S_t"].append(c[2]); out["n_edges"].append(c[3])
            out["mean_tau"].append(c[4] / c[3]); out["bin"].append(c[0]); out["V_aa"].append(c[5]); out["V_t"].append(c[6])
    o = {k: np.asarray(v) for k, v in out.items()}
    o["censored"] = o["S_t"] <= 0
    o["S_t"] = np.where(o["censored"], 0.005, o["S_t"])
    if drop_bin5:
        m = o["bin"] != 4
        o = {k: v[m] for k, v in o.items()}
    o["n_aa_nonpos"] = n_aa_nonpos
    return o


def cell_design(cells, x_center=None, z_center=None):
    y = np.log(cells["S_t"] / cells["S_aa"])
    lx = np.log(cells["mean_tau"])
    lz = np.log(cells["S_aa"] / cells["n_edges"])
    xc = lx.mean() if x_center is None else x_center
    zc = lz.mean() if z_center is None else z_center
    return y, lx - xc, lz - zc, xc, zc


def _try_fit(y, exog, groups, exog_re, start):
    from statsmodels.regression.mixed_linear_model import MixedLM
    from statsmodels.tools.sm_exceptions import ConvergenceWarning
    for meth in ("lbfgs", "bfgs", "powell"):
        for st in ([start, None] if start is not None else [None]):
            with warnings.catch_warnings(record=True) as wl:
                warnings.simplefilter("always")
                try:
                    r = MixedLM(y, exog, groups=groups, exog_re=exog_re).fit(reml=True, method=[meth], start_params=st)
                except Exception:
                    continue
            if bool(r.converged) and np.all(np.isfinite(r.fe_params)) and not any(
                    issubclass(w.category, ConvergenceWarning) for w in wl):
                return r
    return None


def fit_mixed(y, x, z, groups, start=None):
    """MixedLM y ~ 1 + x + z, RE (1 + x | hog), REML, unweighted (trap #6); optimiser chain lbfgs -> bfgs ->
    powell; fallback to a random intercept (flagged) when the slope model does not converge."""
    exog = np.column_stack([np.ones_like(x), x, z])
    flags = []
    r = _try_fit(y, exog, groups, np.column_stack([np.ones_like(x), x]), start)
    conv = r is not None
    if not conv:
        flags.append("re_slope_fallback")
        r = _try_fit(y, exog, groups, None, None)
        if r is None:
            flags.append("fit_failed")
            return {"mu": np.nan, "beta": np.nan, "gamma": np.nan, "flags": flags, "se_model": np.full(3, np.nan),
                    "params": None}
    if conv:
        # REML-maximum polish (pre-confirm fix): statsmodels' optimiser chain accepts the first fit flagged
        # 'converged', which can stop short of the REML maximum (tune dry run: llf -1961.34 vs -1960.62). Restart
        # statsmodels from the fast-REML optimum and keep whichever statsmodels fit has the higher REML llf.
        try:
            from statsmodels.regression.mixed_linear_model import MixedLM, MixedLMParams
            fm = fast_mixed(y, x, z, groups)
            if fm["converged"]:
                Lf = np.zeros((2, 2)); Lf[np.tril_indices(2)] = fm["theta"]
                sp = MixedLMParams.from_components(fe_params=np.array([fm["mu"], fm["beta"], fm["gamma"]]),
                                                   cov_re=Lf @ Lf.T)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    r2 = MixedLM(y, exog, groups=groups, exog_re=np.column_stack([np.ones_like(x), x])).fit(
                        reml=True, method=["lbfgs"], start_params=sp)
                if bool(r2.converged) and np.all(np.isfinite(r2.fe_params)) and r2.llf > r.llf + 1e-6:
                    flags_polish = f"reml_polished:dllf={r2.llf - r.llf:.3g}"
                    r = r2
                    flags.append(flags_polish)
        except Exception:
            pass
    fe = np.asarray(r.fe_params)
    return {"mu": fe[0], "beta": fe[1], "gamma": fe[2], "flags": flags, "se_model": np.asarray(r.bse_fe)[:3],
            "params": r.params if conv else None, "llf": float(r.llf)}


def run_E4(bf, panel: PanelTree, tau_misfit=None, B=200, sensitivities=True, B_sens=None, ages_override=None,
           perm_rng=None, extra_sens=None, input_md5=None, verbose=False, repair_draws=True, engine="fast",
           tau_perm=None, cell_rule="info", info_T=4.0, tip_var=0.0):
    """E4: y = log(S_3Di/S_AA) per HOG x age-bin cell; MixedLM; pigeonhole + age-redraw bootstrap.
    ages_override: (age, lo, hi) arrays per panel node (e.g. permuted ages for V6 / deep-age sensitivities).
    cell_rule / info_T / tip_var: cell merge rule (DEVIATION D30: 'info' = plan AA rule + AA-only predicted-3Di
    information rule; 'aa' = the plan's rule, always reported as sensitivity 'cells_plan_rule'). tip_var = per-protein
    3Di tip-noise variance (tip_sd^2) added to terminal edges in the predicted cell variance."""
    t0 = time.time()
    _bc = build_cells
    _ckw = {"rule": cell_rule, "info_T": info_T, "tip_var": tip_var}

    def _cells(*a, **k):           # every cell construction in this endpoint uses the same rule
        return _bc(*a, **(_ckw | k))
    if ages_override is not None:
        panel = panel.with_ages(*ages_override)
    ev = edge_view(bf, panel, tau_misfit=tau_misfit)
    H = len(ev.hog_ids)
    cells = _cells(ev, panel.age, tau_perm=tau_perm)
    y, x, z, xc, zc = cell_design(cells)
    fit = fit_mixed(y, x, z, cells["hog"])
    flags = list(fit["flags"])
    fchk = fast_mixed(y, x, z, cells["hog"])
    fast_vs_sm = abs(fchk["beta"] - fit["beta"])
    if engine == "fast" and not [f for f in fit["flags"] if not f.startswith("reml_polished")] and fast_vs_sm > 1e-5:
        flags.append(f"fast_reml_mismatch:{fast_vs_sm:.2e}")
    cens_frac = float(cells["censored"].mean()) if len(y) else np.nan
    if cens_frac > 0.02:
        flags.append(f"censored_cells_gt2pct:{cens_frac:.3f}")
    rng = stage_rng("bootstrap_E4", 0)
    pe = panel.panel_edge_ids()
    reps, reps_pooled, rep_flags = [], [], 0
    pooled_pt = pooled_slope(cells)

    def one_rep(evv, drop_bin5=False, internal=False, drop_cens=False):
        W = multinomial_weights(rng, len(evv.hog_ids))
        L = np.zeros(panel.n_nodes); L[pe] = multinomial_weights(rng, len(pe))
        age_b = draw_ages(panel, rng, repair=repair_draws)
        # build cells on the expanded HOG set: edges of HOG h are reused for each copy; Lambda weights inside sums
        ew = _edge_w(evv, np.ones(len(evv.hog_ids)), L)
        c = _cells(evv, age_b, w=ew, drop_bin5=drop_bin5, tau_perm=tau_perm if evv is ev else None)
        if drop_cens:
            m = ~c["censored"]; c = {k: (v[m] if isinstance(v, np.ndarray) else v) for k, v in c.items()}
        # expand cells by W_h copies with new ids
        order = np.argsort(c["hog"], kind="stable")
        ch = c["hog"][order]
        st = np.searchsorted(ch, np.arange(len(evv.hog_ids)), "left"); en = np.searchsorted(ch, np.arange(len(evv.hog_ids)), "right")
        idx, grp = [], []
        gid = 0
        for h in np.flatnonzero(W > 0):
            r_ = order[st[h]:en[h]]
            for _ in range(int(W[h])):
                idx.append(r_); grp.append(np.full(len(r_), gid)); gid += 1
        idx = np.concatenate(idx); grp = np.concatenate(grp)
        cc = {k: v[idx] for k, v in c.items() if isinstance(v, np.ndarray)}
        yy, xx, zz, _, _ = cell_design(cc)
        if engine == "fast":
            f = fast_mixed(yy, xx, zz, grp)
            f["flags"] = [] if f["converged"] else ["re_slope_fallback"]
            if not f["converged"]:
                f = fast_mixed(yy, xx, zz, grp, slope=False); f["flags"] = ["re_slope_fallback"]
        else:
            f = fit_mixed(yy, xx, zz, grp, start=fit.get("params"))
            f["flags"] = [q for q in f["flags"] if not q.startswith("reml_polished")]
        return f, cc

    for b in range(B):
        f, cc = one_rep(ev)
        reps.append(f["beta"]); rep_flags += bool(f["flags"])
        reps_pooled.append(pooled_slope(cc))
        if verbose and b % 20 == 0:
            print(f"   E4 boot {b} {time.time() - t0:.0f}s", flush=True)
    reps = np.asarray(reps)
    n_failed = int((~np.isfinite(reps)).sum())
    reps = reps[np.isfinite(reps)]
    se = float(np.std(reps, ddof=1))
    beta = fit["beta"]
    p = float(2 * (1 - norm.cdf(abs(beta) / se)))
    ci = normal_ci(beta, se); cip = pct_ci(reps)
    if (ci[0] > 0 or ci[1] < 0) != (cip[0] > 0 or cip[1] < 0):
        flags.append("ci_disagree")
    sens = {}
    if sensitivities:
        Bs = B_sens or max(20, B // 5)
        sp = np.asarray(reps_pooled)
        sens["pooled_ratio_of_sums_slope"] = {"estimate": pooled_pt, "ci98": normal_ci(pooled_pt, np.std(sp, ddof=1)), "decision": None}
        sens["wls_1_over_vhg_clustered"] = wls_sensitivity(bf, ev, panel, cells)
        for name, kw in (("censored_dropped", dict(drop_cens=True)), ("bin5_excluded", dict(drop_bin5=True)),
                         ("internal_edges_only", dict(internal=True))):
            evv = edge_view(bf, panel, tau_misfit=tau_misfit, internal_only=True) if kw.get("internal") else ev
            c0 = _cells(evv, panel.age, drop_bin5=kw.get("drop_bin5", False))
            if kw.get("drop_cens"):
                m = ~c0["censored"]; c0 = {k: (v[m] if isinstance(v, np.ndarray) else v) for k, v in c0.items()}
            y0, x0, z0, _, _ = cell_design(c0)
            f0 = fit_mixed(y0, x0, z0, c0["hog"])
            rr = [one_rep(evv, drop_bin5=kw.get("drop_bin5", False), drop_cens=kw.get("drop_cens", False))[0]["beta"] for _ in range(Bs)]
            s0 = float(np.std(rr, ddof=1))
            sens[name] = {"estimate": f0["beta"], "ci98": normal_ci(f0["beta"], s0),
                          "p_unadjusted": float(2 * (1 - norm.cdf(abs(f0["beta"]) / s0))), "B": Bs}
        if cell_rule != "aa":          # the plan's cell rule (known null bias, agents/bias_fix_report.md)
            c0 = _bc(ev, panel.age)
            y0, x0, z0, _, _ = cell_design(c0)
            f0 = fit_mixed(y0, x0, z0, c0["hog"])
            _ckw_save = dict(_ckw); _ckw.update(rule="aa")
            rr = [one_rep(ev)[0]["beta"] for _ in range(Bs)]
            _ckw.clear(); _ckw.update(_ckw_save)
            s0 = float(np.std(rr, ddof=1))
            sens["cells_plan_rule"] = {"estimate": f0["beta"], "ci98": normal_ci(f0["beta"], s0),
                                       "p_unadjusted": float(2 * (1 - norm.cdf(abs(f0["beta"]) / s0))), "B": Bs,
                                       "n_cells": int(len(y0)), "censored_frac": float(c0["censored"].mean()),
                                       "note": "ANALYSIS_PLAN §2.6 rule; biased under the null (D30)"}
    for name, v in (extra_sens or {}).items():
        sens[name] = v
    rec = endpoint_record(
        endpoint="E4", estimand="beta: slope of y=log(S_3Di/S_AA) on centred log mean cell age (MixedLM); beta_k = -beta",
        estimate=beta, se_boot=se, ci98=ci, ci_pct98=cip,
        one_sided_bound={"lower99": beta - Z98 * se, "upper99": beta + Z98 * se},
        p_value=p, p_holm=None, decision="pending_holm", n_hogs=H, n_units=int(len(y)),
        n_excluded=ev.n_excluded | {"cells_aa_nonpositive": int(cells["n_aa_nonpos"])}, sensitivity=sens,
        alignment_robust=None, flags=flags, input_md5=input_md5 or md5_arrays(ev.b_aa, ev.b_t),
        beta_k=-beta, mu=fit["mu"], gamma=fit["gamma"], x_center=xc, z_center=zc, se_model=fit["se_model"],
        censored_frac=cens_frac, cell_rule=cell_rule, info_T=info_T, tip_var=tip_var, B=B, rep_fallback_count=rep_flags, rep_failed=n_failed, mixed_engine_boot=engine, fast_vs_statsmodels_abs_diff=fast_vs_sm, runtime_s=time.time() - t0)
    return rec


def pooled_slope(cells):
    """Pooled per-bin ratio-of-sums: log(sum S_t / sum S_aa) per bin on centred log mean tau (OLS)."""
    b = cells["bin"]
    ub = np.unique(b)
    if len(ub) < 2:
        return np.nan
    R = np.array([np.log(cells["S_t"][b == u].sum() / cells["S_aa"][b == u].sum()) for u in ub])
    T = np.array([np.log(np.average(cells["mean_tau"][b == u], weights=cells["n_edges"][b == u])) for u in ub])
    T = T - T.mean()
    return float((T * (R - R.mean())).sum() / (T * T).sum())


def wls_sensitivity(bf, ev, panel, cells):
    """WLS with weights 1/v_hg (delta method from the sandwich scores) and HOG-clustered SEs."""
    import statsmodels.api as sm
    tau = 0.5 * (panel.age[ev.child] + panel.age[ev.parent])
    g = np.clip(np.searchsorted(BIN_EDGES, tau, side="right") - 1, 0, 4)
    hf_by = {hf.hog: hf for hf in bf.hogs}
    ka = bf.metrics.index("aa_ml"); kt = bf.metrics.index("tdi_ml")
    v = np.full(len(cells["hog"]), np.nan)
    # merged cells: approximate cell membership by the cell's own bin and all bins merged into it is not stored;
    # use the bins of edges whose tau-bin equals the cell bin (conservative approximation, documented)
    for i, (h, bn) in enumerate(zip(cells["hog"], cells["bin"])):
        hf = hf_by[ev.hog_ids[h]]
        m = (ev.hog_code == h) & (g == bn)
        ht = hf.ht
        a = np.zeros(ht.n_params)
        # edges of this HOG in view order == ht edge order (internal_only views excluded here)
        sel = np.flatnonzero(ev.hog_code == h)
        mm = m[sel]
        if len(mm) != len(ht.param_of_edge):
            continue
        np.add.at(a, ht.param_of_edge[mm], ht.share[mm])
        ga = a @ hf.fit["G"][ka]; gt = a @ hf.fit["G"][kt]
        SA, ST = cells["S_aa"][i], cells["S_t"][i]
        v[i] = (gt @ gt) / ST ** 2 + (ga @ ga) / SA ** 2 - 2 * (ga @ gt) / (SA * ST)
    y, x, z, _, _ = cell_design(cells)
    ok = np.isfinite(v) & (v > 0)
    X = sm.add_constant(np.column_stack([x, z]))
    r = sm.WLS(y[ok], X[ok], weights=1 / v[ok]).fit(cov_type="cluster", cov_kwds={"groups": cells["hog"][ok]})
    b, s = float(r.params[1]), float(r.bse[1])
    return {"estimate": b, "ci98": normal_ci(b, s), "p_unadjusted": float(r.pvalues[1]), "n_cells": int(ok.sum())}
