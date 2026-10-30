"""E5 (ANALYSIS_PLAN §4.5): lineage ratios, parametric-bootstrap floor F, excess lineage variance sigma2_L,
disattenuated sister-contrast slope beta, spike-in power, N0 decision with guards."""
from __future__ import annotations

import time

import numpy as np
import scipy.sparse as sps
from scipy.stats import norm

from .common import Z98, endpoint_record, normal_ci, pct_ci, stage_rng, md5_arrays
from .bootstrap import multinomial_weights
from .endpoints_branch import edge_view, EdgeView

CENSOR = 0.005


class E5Data:
    """Precomputed structures for E5 on one EdgeView."""

    def __init__(self, bf, panel, ev: EdgeView, fe: dict, species_plddt=None, min_hogs=50, min_D=0.05,
                 min_hogs_contrast=50, internal_only=False, min_t_snr=None):
        self.bf, self.panel, self.ev, self.fe = bf, panel, ev, fe
        H = len(ev.hog_ids)
        self.H = H
        # ---- lineages (simple HOG edges mapped to panel edge e)
        s = ev.simple & (~ev.is_terminal if internal_only else True)
        rows = np.flatnonzero(s)
        self.internal_only = internal_only
        child = ev.child[rows]
        hc = ev.hog_code[rows]
        uniq = np.unique(child)
        nh = np.array([len(np.unique(hc[child == c])) for c in uniq])
        D = np.array([ev.b_aa[rows][child == c].sum() for c in uniq])
        keep_l = (nh >= min_hogs) & (D >= min_D)
        self.min_t_snr = min_t_snr
        if min_t_snr is not None:
            # PROPOSED DEVIATION D8: also require expected 3Di-sum precision r_bar*D_e / sqrt(sum SE_3Di^2) >= min_t_snr
            # (selection on AA sums and SEs only, never on observed 3Di values)
            rbar = ev.b_t.sum() / ev.b_aa.sum()
            se_t = np.array([np.sqrt(np.sum(ev.se_t[rows][child == c] ** 2)) for c in uniq])
            keep_l &= (rbar * D / np.maximum(se_t, 1e-12)) >= min_t_snr
        lin = uniq[keep_l]
        self.lineages = lin
        pos = {int(c): i for i, c in enumerate(lin)}
        m = np.isin(child, lin)
        self.l_rows = rows[m]
        self.l_idx = np.array([pos[int(c)] for c in child[m]], dtype=int)
        self.l_hog = hc[m]
        nL = len(lin)
        age = panel.age
        tau_e = 0.5 * (age[lin] + age[panel.parent[lin]])
        self.x_e = np.log(tau_e) - fe["x_center"]
        # pLDDT of descendant species
        if species_plddt is not None and nL:
            pl = []
            for c in lin:
                leaves = [panel.taxid[panel.leaf_order[k]] for k in range(len(panel.leaf_order)) if (panel.clade[c] >> k) & 1]
                v = [species_plddt.get(int(t), np.nan) for t in leaves]
                pl.append(np.nanmean(v) if np.isfinite(v).any() else np.nan)   # panel species outside the sample: ignored
            self.plddt = np.asarray(pl, float)
            if not np.isfinite(self.plddt).all():                       # no sampled descendant: mean-impute (driver flags)
                self.plddt = np.where(np.isfinite(self.plddt), self.plddt, np.nanmean(self.plddt) if np.isfinite(self.plddt).any() else 0.0)
        else:
            self.plddt = None
        # ---- sister contrasts: sparse path-length operators
        self._build_contrasts(min_hogs_contrast)

    # ------------------------------------------------------------------ lineage statistic
    def sigma2_obs(self, b_aa, b_t, W=None, lam=None):
        """Return (sigma2_obs, ell, logR, n_censored) for edge values b_aa, b_t (arrays over ev rows or (P, rows))."""
        nL = len(self.lineages)
        W = np.ones(self.H) if W is None else W
        lam = np.ones(nL) if lam is None else lam
        w = W[self.l_hog]
        den = np.bincount(self.l_idx, w * b_aa[self.l_rows], nL)
        num = np.bincount(self.l_idx, w * b_t[self.l_rows], nL)
        cnt = np.bincount(self.l_idx, w, nL)
        cens = num <= 0
        num = np.where(cens, CENSOR, num)
        ok = (den > 0) & (cnt > 0) & (lam > 0)
        logR = np.full(nL, np.nan)
        logR[ok] = np.log(num[ok] / den[ok])
        z_e = np.full(nL, np.nan)
        z_e[ok] = np.log(den[ok] / cnt[ok]) - self.fe["z_center"]
        ell = logR - (self.fe["mu"] + self.fe["beta"] * self.x_e + self.fe["gamma"] * z_e)
        lw = lam[ok]
        e = ell[ok]
        if self.plddt is not None:
            X = np.column_stack([np.ones(ok.sum()), self.plddt[ok]])
            sw = np.sqrt(lw)
            coef = np.linalg.lstsq(X * sw[:, None], e * sw, rcond=None)[0]
            rho = e - X @ coef
        else:
            rho = e - np.average(e, weights=lw)
        mu = np.average(rho, weights=lw)
        s2 = np.sum(lw * (rho - mu) ** 2) / max(lw.sum() - 1, 1)
        return s2, ell, logR, int(cens.sum())

    # ------------------------------------------------------------------ contrasts
    def _build_contrasts(self, min_hogs):
        panel, ev = self.panel, self.ev
        hf_by = {hf.hog: hf for hf in self.bf.hogs}
        # ev rows per HOG
        order = np.argsort(ev.hog_code, kind="stable")
        starts = np.searchsorted(ev.hog_code[order], np.arange(self.H), "left")
        ends = np.searchsorted(ev.hog_code[order], np.arange(self.H), "right")
        contrasts = []
        for v in np.flatnonzero(~panel.is_leaf):
            kids = sorted(panel.children[v], key=lambda c: panel.label[c])
            for a, b in zip(kids[:-1], kids[1:]):
                contrasts.append((v, a, b))
        rec_rows, rec_cols, ca, cb = [], [], [], []
        rec_meta = []
        nrec = 0
        leafset = lambda node: panel.clade[node]
        for h in range(self.H):
            hf = hf_by[ev.hog_ids[h]]
            ht = hf.ht
            evrows = order[starts[h]:ends[h]]
            if len(evrows) != len(ht.edge_child):
                continue   # internal-only views: contrasts are computed on full views only
            sp_bits = np.array([1 << panel.leaf_pos[panel.leaf_of_taxid[int(t)]] for t in ht.species], dtype=object)
            mem = ht.member
            sizes = mem.sum(0)
            for ci, (v, a, b) in enumerate(contrasts):
                inA = np.array([(int(sb) & leafset(a)) != 0 for sb in sp_bits])
                inB = np.array([(int(sb) & leafset(b)) != 0 for sb in sp_bits])
                if not inA.any() or not inB.any():
                    continue
                need = inA | inB
                cand = np.flatnonzero(mem[need].all(0))
                if cand.size:
                    u = cand[np.argmin(sizes[cand])]
                    umem = mem[:, u]
                else:
                    u = -1
                    if ht.root_merged:
                        continue      # split of the HOG root edge is not identifiable (D7)
                    umem = np.ones(len(ht.species), bool)
                below = np.all(mem <= umem[:, None], axis=0) & (mem.sum(0) < umem.sum())
                if self.internal_only:
                    below = below & ~ht.is_terminal
                cA = (mem[inA][:, below]).sum(0) / inA.sum()
                cB = (mem[inB][:, below]).sum(0) / inB.sum()
                cols = evrows[np.flatnonzero(below)]
                rec_rows += [nrec] * len(cols); rec_cols += list(cols); ca += list(cA); cb += list(cB)
                rec_meta.append((ci, h))
                nrec += 1
        nE = len(ev.child)
        self.contrasts = contrasts
        meta = np.array(rec_meta, dtype=int).reshape(-1, 2)
        self.c_rec_contrast = meta[:, 0]
        self.c_rec_hog = meta[:, 1]
        self.CA = sps.csr_matrix((ca, (rec_rows, rec_cols)), shape=(nrec, nE))
        self.CB = sps.csr_matrix((cb, (rec_rows, rec_cols)), shape=(nrec, nE))
        nh = np.bincount(self.c_rec_contrast, minlength=len(contrasts))
        keep_c = np.flatnonzero(nh >= min_hogs)
        self.c_keep = keep_c
        self.c_map = -np.ones(len(contrasts), int)
        self.c_map[keep_c] = np.arange(len(keep_c))
        self.c_rec_k = self.c_map[self.c_rec_contrast]

    def contrast_xy(self, b_aa, b_t, W=None):
        """x_v, y_v for kept contrasts."""
        W = np.ones(self.H) if W is None else W
        nk = len(self.c_keep)
        m = self.c_rec_k >= 0
        w = W[self.c_rec_hog][m]
        k = self.c_rec_k[m]
        LAa = (self.CA @ b_aa)[m]; LBa = (self.CB @ b_aa)[m]
        LAt = (self.CA @ b_t)[m]; LBt = (self.CB @ b_t)[m]
        sA = np.bincount(k, w * LAa, nk); sB = np.bincount(k, w * LBa, nk)
        tA = np.bincount(k, w * LAt, nk); tB = np.bincount(k, w * LBt, nk)
        with np.errstate(invalid="ignore", divide="ignore"):
            x = np.log(sA / sB); y = np.log(tA / tB)
        return x, y


def _beta(x, y, c, vx, wv=None):
    wv = np.ones(len(x)) if wv is None else wv
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(c) & np.isfinite(vx)
    return float((np.sum(wv[ok] * x[ok] * y[ok]) - np.sum(wv[ok] * c[ok])) /
                 (np.sum(wv[ok] * x[ok] ** 2) - np.sum(wv[ok] * vx[ok])))


def parametric_floor(dat: E5Data, s2=None, rho=0.0, P=100, nu_t=0.0, rng=None):
    """F = mean sigma2_obs over P parametric replicates simulated under N0 from the fitted tree
    (b_3Di = r_h exp(beta (x_e - xbar_h)) b_AA), bivariate noise Var = SE_cal^2 + s2, corr rho; branch refit is
    exact-linear (b* = A d*) with the observed second-pass weights."""
    rng = np.random.default_rng(rng)
    ev, panel, bf = dat.ev, dat.panel, dat.bf
    hf_by = {hf.hog: hf for hf in bf.hogs}
    ka = bf.metrics.index("aa_ml"); kt = bf.metrics.index("tdi_ml")
    order = np.argsort(ev.hog_code, kind="stable")
    starts = np.searchsorted(ev.hog_code[order], np.arange(dat.H), "left")
    ends = np.searchsorted(ev.hog_code[order], np.arange(dat.H), "right")
    nE = len(ev.child)
    BA = np.zeros((P, nE)); BT = np.zeros((P, nE))
    beta = dat.fe["beta"]
    for h in range(dat.H):
        hf = hf_by[ev.hog_ids[h]]
        rows = order[starts[h]:ends[h]]
        ht = hf.ht
        bA = ev.b_aa[rows]
        r_h = ev.b_t[rows].sum() / bA.sum() if bA.sum() > 0 else np.nan
        if not np.isfinite(r_h) or r_h <= 0:
            r_h = np.nanmedian(ev.b_t) / np.nanmedian(ev.b_aa) if np.nanmedian(ev.b_aa) > 0 else 0.3
        tau = 0.5 * (panel.age[ht.edge_child] + panel.age[ht.edge_parent])
        xe = np.log(np.maximum(tau, 1e-3))
        bT = r_h * np.exp(beta * (xe - xe.mean())) * bA
        pa = np.zeros(ht.n_params); pt = np.zeros(ht.n_params)
        np.add.at(pa, ht.param_of_edge, bA); np.add.at(pt, ht.param_of_edge, bT)
        dA = hf.X @ pa
        dT = hf.X @ pt + 2 * nu_t
        okA, okT = hf.ok[ka], hf.ok[kt]
        s2a = hf.fit["s2"][ka] if s2 is None else s2["aa_ml"]
        s2t = hf.fit["s2"][kt] if s2 is None else s2["tdi_ml"]
        va = hf.se[ka] ** 2 + s2a
        vt = hf.se[kt] ** 2 + s2t
        z1 = rng.standard_normal((P, len(dA)))
        z2 = rho * z1 + np.sqrt(1 - rho * rho) * rng.standard_normal((P, len(dA)))
        dsA = dA + np.sqrt(np.where(okA, va, 0)) * z1
        dsT = dT + np.sqrt(np.where(okT, vt, 0)) * z2
        pA = dsA @ hf.fit["A"][ka].T
        pT = dsT @ hf.fit["A"][kt].T
        BA[:, rows] = ht.share * pA[:, ht.param_of_edge]
        BT[:, rows] = ht.share * pT[:, ht.param_of_edge] - np.where(ht.is_terminal, nu_t, 0.0)
    vals = np.array([dat.sigma2_obs(BA[p], BT[p])[0] for p in range(P)])
    return float(np.mean(vals)), vals


def e5_core(dat: E5Data, b_aa, b_t, F, B=1000, B_cv=200, rng=None, cv=None):
    """sigma2_L, beta_contrast with bootstrap CIs (F fixed; c_v, vx_v from B_cv HOG-pigeonhole replicates
    unless given)."""
    rng = np.random.default_rng(rng)
    s2o, ell, logR, ncens = dat.sigma2_obs(b_aa, b_t)
    s2L = s2o - F
    x, y = dat.contrast_xy(b_aa, b_t)
    if cv is None:
        xs, ys = [], []
        for _ in range(B_cv):
            W = multinomial_weights(rng, dat.H)
            xb, yb = dat.contrast_xy(b_aa, b_t, W)
            xs.append(xb); ys.append(yb)
        xs, ys = np.array(xs), np.array(ys)
        c_v = np.nanmean((xs - np.nanmean(xs, 0)) * (ys - np.nanmean(ys, 0)), 0) * B_cv / max(B_cv - 1, 1)
        vx_v = np.nanvar(xs, 0, ddof=1)
        cv = (c_v, vx_v)
    c_v, vx_v = cv
    beta = _beta(x, y, c_v, vx_v) if len(x) else np.nan
    nL = len(dat.lineages)
    s2b, bb, bbh = [], [], []
    for _ in range(B):
        W = multinomial_weights(rng, dat.H)
        lam = multinomial_weights(rng, nL) if nL else None
        s2b.append(dat.sigma2_obs(b_aa, b_t, W, lam)[0] - F)
        if len(x):
            xb, yb = dat.contrast_xy(b_aa, b_t, W)
            wv = multinomial_weights(rng, len(xb))     # contrast-level pigeonhole weights (clarification C9)
            bb.append(_beta(xb, yb, c_v, vx_v, wv))
            bbh.append(_beta(xb, yb, c_v, vx_v))
    s2b = np.array(s2b); bb = np.array(bb); bbh = np.array(bbh)
    se_s = float(np.nanstd(s2b, ddof=1)); se_b = float(np.nanstd(bb, ddof=1)) if len(bb) else np.nan
    se_bh = float(np.nanstd(bbh, ddof=1)) if len(bbh) else np.nan
    return {"sigma2_obs": s2o, "sigma2_L": s2L, "sigma2_L_se": se_s,
            "sigma2_L_lower99": s2L - Z98 * se_s, "sigma2_L_upper99": s2L + Z98 * se_s,
            "sigma_L": float(np.sqrt(max(0.0, s2L))), "sigma_L_upper99": float(np.sqrt(max(0.0, s2L + Z98 * se_s))),
            "sigma_L_lower99": float(np.sqrt(max(0.0, s2L - Z98 * se_s))),
            "beta": beta, "beta_se": se_b, "beta_ci98": normal_ci(beta, se_b), "beta_ci_pct98": pct_ci(bb),
            "beta_se_hog_only": se_bh, "beta_ci98_hog_only": normal_ci(beta, se_bh),
            "sigma2_L_ci_pct98": pct_ci(s2b), "n_lineages": nL, "n_contrasts": int(np.sum(np.isfinite(x))),
            "n_censored_lineages": ncens, "cv": cv, "ell": ell}


def e5_margins(core):
    ci = core["beta_ci98"]
    beta_in = ci[0] >= 0.8 and ci[1] <= 1.25
    beta_out = ci[1] < 0.8 or ci[0] > 1.25
    sig_ok = core["sigma2_L_upper99"] < 0.01
    sig_ex = core["sigma2_L_lower99"] > 0.01
    return beta_in, beta_out, sig_ok, sig_ex


def e5_pvalues(core):
    b, sb = core["beta"], core["beta_se"]
    s, ss = core["sigma2_L"], core["sigma2_L_se"]
    p_b_lt = 1 - norm.cdf((b - 0.8) / sb)          # H0: beta <= 0.8
    p_b_gt = norm.cdf((b - 1.25) / sb)             # H0: beta >= 1.25
    p_s = norm.cdf((s - 0.01) / ss)                # H0: sigma2_L >= 0.01
    rej = {"beta_lt_0.8": norm.cdf((b - 0.8) / sb), "beta_gt_1.25": 1 - norm.cdf((b - 1.25) / sb),
           "sigma_gt_0.10": 1 - norm.cdf((s - 0.01) / ss)}
    return float(max(p_b_lt, p_b_gt, p_s)), {k: float(v) for k, v in rej.items()}


def e5_decide(core, guards=None, e6_significant=None, power=None):
    beta_in, beta_out, sig_ok, sig_ex = e5_margins(core)
    p_iut, p_rej = e5_pvalues(core)
    if beta_out or sig_ex:
        guard_ok = all((e5_margins(g)[1] or e5_margins(g)[3]) for g in (guards or {}).values())
        dec = "N0_rejected" if guard_ok else "artefact_suspect"
        rej_ps = [p_rej["beta_lt_0.8"] if core["beta_ci98"][1] < 0.8 else None,
                  p_rej["beta_gt_1.25"] if core["beta_ci98"][0] > 1.25 else None,
                  p_rej["sigma_gt_0.10"] if sig_ex else None]
        p = min(v for v in rej_ps if v is not None)
    elif beta_in and sig_ok and bool(e6_significant) and (power is not None and power >= 0.80):
        dec, p = "N0_accepted", p_iut
    else:
        dec, p = "inconclusive", p_iut
    return dec, float(p), {"margins_met": bool(beta_in and sig_ok), "beta_in": beta_in, "beta_out": beta_out,
                           "sigma_ok": sig_ok, "sigma_excess": sig_ex, "p_iut": p_iut, "p_reject": p_rej,
                           "e6_significant": e6_significant, "power": power}


def spikein_power(dat: E5Data, F, cv, f=1.5, reps=50, B=200, frac=0.10, stage_extra=0):
    """Spike-in power (ANALYSIS_PLAN §4.5): multiply b_3Di on simple edges of a random 10% of lineages by f
    (exactly what adding (f-1) b_3Di,he to every crossing pair and refitting gives, since A X = I), run E5,
    detection = sigma2_L lower 99% > 0 or beta 98% CI excludes 1. Returns (power, accept_rate, details)."""
    ev = dat.ev
    nL = len(dat.lineages)
    det, acc = [], []
    for r in range(reps):
        rng = stage_rng("spikein", stage_extra, r)
        chosen = rng.choice(nL, size=max(1, int(round(frac * nL))), replace=False)
        bt = ev.b_t.copy()
        sel = dat.l_rows[np.isin(dat.l_idx, chosen)]
        bt[sel] *= f
        core = e5_core(dat, ev.b_aa, bt, F, B=B, rng=rng, cv=cv)
        d = (core["sigma2_L_lower99"] > 0) or not (core["beta_ci98"][0] <= 1 <= core["beta_ci98"][1])
        bi, bo, so, se_ = e5_margins(core)
        det.append(bool(d)); acc.append(bool(bi and so))
    return float(np.mean(det)), float(np.mean(acc)), det


def run_E5(bf, panel, e4_rec, tau_misfit=None, species_plddt=None, rho=0.0, nu_t=0.0, P=100, B=1000, B_cv=200,
           spike_reps=50, spike_B=200, spike_fs=(1.5, 1.25, 2.0, 1.0), guards=None, e6_significant=None,
           min_hogs=50, input_md5=None, power_override=None, min_t_snr=None):
    """Full E5. `guards` = {name: BranchFit} variants (e.g. hidden paralogy removed, single-copy HOGs); the
    internal-edges-only guard is computed here."""
    t0 = time.time()
    fe = {k: e4_rec[k] for k in ("mu", "gamma", "x_center", "z_center")}
    fe["beta"] = e4_rec["estimate"]
    ev = edge_view(bf, panel, tau_misfit=tau_misfit)
    dat = E5Data(bf, panel, ev, fe, species_plddt, min_hogs=min_hogs, min_hogs_contrast=min_hogs, min_t_snr=min_t_snr)
    F, Fvals = parametric_floor(dat, rho=rho, P=P, nu_t=nu_t, rng=stage_rng("bootstrap_E5", 1))
    core = e5_core(dat, ev.b_aa, ev.b_t, F, B=B, B_cv=B_cv, rng=stage_rng("bootstrap_E5", 0))
    power = {}
    for i, f in enumerate(spike_fs):
        pw, ac, _ = spikein_power(dat, F, core["cv"], f=f, reps=spike_reps, B=spike_B, stage_extra=i)
        power[f"f={f}"] = {"detect_rate": pw, "accept_rate": ac}
    pw15 = power.get("f=1.5", {}).get("detect_rate") if power_override is None else power_override
    gres = {}
    # internal-edges-only guard
    dati = E5Data(bf, panel, ev, fe, species_plddt, min_hogs=min_hogs, min_hogs_contrast=min_hogs, internal_only=True,
                  min_t_snr=min_t_snr)
    if len(dati.lineages) >= 3:
        Fi, _ = parametric_floor(dati, rho=rho, P=max(20, P // 5), nu_t=nu_t, rng=stage_rng("bootstrap_E5", 2))
        gres["internal_edges_only"] = e5_core(dati, ev.b_aa, ev.b_t, Fi, B=max(100, B // 5),
                                              rng=stage_rng("bootstrap_E5", 3))
    for name, bfg in (guards or {}).items():
        evg = edge_view(bfg, panel, tau_misfit=tau_misfit)
        dg = E5Data(bfg, panel, evg, fe, species_plddt, min_hogs=min_hogs, min_hogs_contrast=min_hogs, min_t_snr=min_t_snr)
        Fg, _ = parametric_floor(dg, rho=rho, P=max(20, P // 5), nu_t=nu_t, rng=stage_rng("bootstrap_E5", 4))
        gres[name] = e5_core(dg, evg.b_aa, evg.b_t, Fg, B=max(100, B // 5), rng=stage_rng("bootstrap_E5", 5))
    dec, p, info = e5_decide(core, gres, e6_significant, pw15)
    flags = []
    if power_override is None and pw15 is not None and pw15 < 0.80:
        flags.append("power_lt_0.80")
    if core["n_censored_lineages"]:
        flags.append(f"censored_lineages:{core['n_censored_lineages']}")
    sens = {k: {"estimate": {"beta": g["beta"], "sigma2_L": g["sigma2_L"]},
                "ci98": {"beta": g["beta_ci98"], "sigma2_L": [g["sigma2_L_lower99"], g["sigma2_L_upper99"]]},
                "decision": "excess" if (e5_margins(g)[1] or e5_margins(g)[3]) else "no_excess"} for k, g in gres.items()}
    rec = endpoint_record(
        endpoint="E5", estimand="N0 (proportional lineage effects): disattenuated sister-contrast slope beta and "
                                "excess lineage SD sigma_L",
        estimate=core["beta"], se_boot=core["beta_se"], ci98=core["beta_ci98"], ci_pct98=core["beta_ci_pct98"],
        one_sided_bound={"sigma2_L_lower99": core["sigma2_L_lower99"], "sigma2_L_upper99": core["sigma2_L_upper99"],
                         "sigma_L_upper99": core["sigma_L_upper99"], "sigma_L_lower99": core["sigma_L_lower99"]},
        p_value=p, p_holm=None, decision=dec, n_hogs=dat.H, n_units=core["n_lineages"],
        n_excluded=ev.n_excluded, sensitivity=sens, alignment_robust=None, flags=flags,
        input_md5=input_md5 or md5_arrays(ev.b_aa, ev.b_t),
        sigma_L=core["sigma_L"], sigma2_L=core["sigma2_L"], sigma2_obs=core["sigma2_obs"], beta_contrast=core["beta"],
        power=power, floor_F=F, floor_F_sd=float(np.std(Fvals, ddof=1)), n_lineages=core["n_lineages"],
        n_contrasts=core["n_contrasts"], decision_detail=info, min_t_snr=min_t_snr,
        beta_ci98_hog_only=core["beta_ci98_hog_only"], B=B, P=P, runtime_s=time.time() - t0)
    return rec, dat, core
