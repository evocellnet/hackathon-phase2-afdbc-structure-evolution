"""Branch engine (ANALYSIS_PLAN §2.5; DESIGN §11.1; R2-M1/M2/M8; R3-M4/M12/M13).

Per HOG h, metric m, alignment: unconstrained two-pass WLS of observed pairwise distances on the path-
incidence design of T_h (panel tree pruned to S_h, unary nodes suppressed), HC1 sandwich SEs, cross-metric
covariance with AA, tip correction nu_m, RSS_rel, NJ-misfit ratio (AA), hidden-paralogy flags.
No clipping anywhere (trap #1); NNLS only for display.

Vectorisation: the design X is shared by all metrics (and alignments, if passed together) of one HOG, so all
K = metrics x alignments fits of a HOG are solved in one batched call; HOGs are looped in python
(~2-4 ms per HOG incl. pruning, NJ and NNLS).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import nnls

from .trees import PanelTree, HogTree, prune, design_matrix, nj_splits, edge_ages

METRIC_COLS = {  # metric -> (value, se, saturation flag)
    "aa_ml": ("aa_ml", "aa_ml_se", "aa_sat"),
    "tdi_ml": ("tdi_ml", "tdi_ml_se", "tdi_sat"),
    "lddt_D": ("lddt_D", "lddt_D_se", "lddt_sat"),
}


# ------------------------------------------------------------------------------------------ core WLS
def wls_two_pass(X: np.ndarray, d: np.ndarray, se: np.ndarray, ok: np.ndarray, min_extra: int = 5,
                 se_method: str = "hc1"):
    """Batched two-pass WLS for K response vectors sharing the design X.

    X (n,p); d, se, ok (K,n). Failed / saturated pairs have ok=False -> weight 0 (they still count in n_pairs).
    Pass 1: w = 1/SE^2 ; s2 = max(0, RSS_w/(n_ok-p) * mean(SE^2) - mean(SE^2)).
    Pass 2: w = 1/(SE^2 + s2) ; b = (X'WX)^-1 X'W d (no clipping).
    HC1 sandwich: V = A diag(r^2) A' * n_ok/(n_ok-p), A = (X'WX)^-1 X'W  (se_method='hc1', the plan).
    se_method='hc3' uses r_i/(1-h_ii); 'model' uses (X'WX)^-1 (= A diag(1/w) A'). Proposed deviation D2.
    The 'scores' G = A diag(r~) are returned so that any linear combination / cross-metric covariance
    can be formed as G1 G2'.
    Returns dict of arrays; fits with n_ok < p + min_extra are flagged underdetermined (b = NaN).
    """
    K, n = d.shape
    p = X.shape[1]
    ok = ok & np.isfinite(d) & np.isfinite(se) & (se > 0)
    n_ok = ok.sum(1)
    under = n_ok < p + min_extra
    # identifiability (pre-confirm fix C/D, presconfirm_fixes_report): the count rule alone admits designs whose
    # ok rows do not identify every edge (e.g. after hidden-paralogy pair removal); X'WX is then singular and b
    # explodes (|b| ~ 1e12 seen on the tune split). Such fits are flagged rank-deficient and treated as underdetermined.
    rank_def = np.zeros(K, bool)
    for k in np.flatnonzero(~under):
        if np.linalg.matrix_rank(X[ok[k]]) < p:
            rank_def[k] = True
    under = under | rank_def
    dz = np.where(ok, d, 0.0)
    se2 = np.where(ok, se * se, 1.0)
    out = {"b": np.full((K, p), np.nan), "A": np.zeros((K, p, n)), "r": np.zeros((K, n)), "w": np.zeros((K, n)),
           "s2": np.full(K, np.nan), "rss_rel": np.full(K, np.nan), "n_ok": n_ok, "under": under.copy(),
           "hc1": np.full(K, np.nan), "V": np.full((K, p, p), np.nan), "G": np.zeros((K, p, n)),
           "se_method": se_method, "rank_def": rank_def}
    act = ~under
    if not act.any():
        return out
    Xa = X[None]

    def solve(w):
        XtW = np.swapaxes(Xa * w[:, :, None], 1, 2)         # (K,p,n)
        XtWX = XtW @ X                                      # (K,p,p)
        try:
            A = np.linalg.solve(XtWX, XtW)                  # (K,p,n)
        except np.linalg.LinAlgError:
            A = np.stack([np.linalg.lstsq(XtWX[k], XtW[k], rcond=None)[0] for k in range(len(w))])
        return A

    w1 = np.where(ok, 1.0 / se2, 0.0)
    w1[~act] = 0.0
    A1 = solve(np.where(act[:, None], w1, 1e-12))
    b1 = np.einsum("kpn,kn->kp", A1, dz)
    r1 = np.where(ok, dz - b1 @ X.T, 0.0)
    rssw = (w1 * r1 * r1).sum(1)
    mse2 = np.array([se2[k, ok[k]].mean() if ok[k].any() else np.nan for k in range(K)])
    with np.errstate(invalid="ignore", divide="ignore"):
        s2 = np.maximum(0.0, rssw / np.maximum(n_ok - p, 1) * mse2 - mse2)
    w2 = np.where(ok, 1.0 / (se2 + s2[:, None]), 0.0)
    A = solve(np.where(act[:, None], w2, 1e-12))
    b = np.einsum("kpn,kn->kp", A, dz)
    r = np.where(ok, dz - b @ X.T, 0.0)
    A = np.where(ok[:, None, :], A, 0.0)
    hc1 = n_ok / np.maximum(n_ok - p, 1)
    lev = np.einsum("np,kpn->kn", X, A)                     # leverages h_ii (0 for weight-0 rows)
    if se_method == "hc1":
        rs = r * np.sqrt(hc1)[:, None]
    elif se_method == "hc3":
        rs = np.where(lev < 1 - 1e-8, r / np.maximum(1 - lev, 1e-8), 0.0)
    elif se_method != "model":
        raise ValueError(se_method)
    if se_method == "model":
        # model-based (X'WX)^-1 with w = 1/(SE^2+s2); scores use the model SD instead of the residual
        rs = np.where(ok, np.sqrt(1.0 / np.maximum(w2, 1e-300)), 0.0)
    G = A * rs[:, None, :]
    V = G @ np.swapaxes(G, 1, 2)
    with np.errstate(invalid="ignore", divide="ignore"):
        rss_rel = (w2 * r * r).sum(1) / (w2 * dz * dz).sum(1)
    for k in np.flatnonzero(act):
        out["b"][k] = b[k]; out["A"][k] = A[k]; out["r"][k] = r[k]; out["w"][k] = w2[k]
        out["s2"][k] = s2[k]; out["rss_rel"][k] = rss_rel[k]; out["hc1"][k] = hc1[k]; out["V"][k] = V[k]
        out["G"][k] = G[k]
    return out


def nnls_fit(X, d, w, ok):
    sw = np.sqrt(np.where(ok, w, 0.0))
    try:
        return nnls(X * sw[:, None], np.where(ok, d, 0.0) * sw, maxiter=50 * X.shape[1])[0]
    except Exception:
        return np.full(X.shape[1], np.nan)


# ------------------------------------------------------------------------------------------ containers
@dataclass
class HogFit:
    hog: str
    split: str
    ht: HogTree
    ia: np.ndarray
    ib: np.ndarray
    X: np.ndarray
    metrics: list
    d: np.ndarray        # (K,n) observed
    se: np.ndarray       # (K,n) calibrated SE
    ok: np.ndarray       # (K,n)
    fit: dict            # output of wls_two_pass
    excluded: dict = field(default_factory=dict)   # metric -> reason ('' if kept)
    misfit_ratio: float = np.nan
    hidden_para: np.ndarray | None = None
    pair_index: np.ndarray | None = None           # row indices into the input pair table


@dataclass
class BranchFit:
    aln: str
    metrics: list
    nu: dict
    hogs: list                      # list[HogFit]
    edges: dict                     # columnar per HOG-edge x metric table (numpy arrays)
    fits: dict                      # columnar per HOG x metric table
    pair_flags: dict                # columnar per pair: hidden paralogy

    def edge_frame(self, panel: PanelTree, with_ages=True):
        import polars as pl
        e = dict(self.edges)
        cols = {k: v for k, v in e.items() if not k.startswith("_")}
        cols["panel_edge"] = pl.Series("panel_edge", [None if v is None else str(v) for v in cols["panel_edge"]], dtype=pl.Utf8)
        for k in ("hog", "split", "metric", "edge_id", "panel_child", "panel_parent"):
            cols[k] = cols[k].astype(str)
        if not with_ages:
            for k in ("t_child", "t_parent", "tau"):
                cols.pop(k, None)
        return pl.DataFrame(cols)

    def fit_frame(self):
        import polars as pl
        cols = {k: v for k, v in self.fits.items() if not k.startswith("_")}
        for k in ("hog", "split", "metric", "reason"):
            cols[k] = cols[k].astype(str)
        return pl.DataFrame(cols)


# ------------------------------------------------------------------------------------------ driver
def _col(df, name):
    if hasattr(df, "get_column"):
        return df.get_column(name).to_numpy()
    return np.asarray(df[name])


def fit_branches(pairs, panel: PanelTree, metrics=("aa_ml", "tdi_ml", "lddt_D"), nu=None, kappa=None,
                 aln="tmalign", hog_species: dict | None = None, tau_misfit: float | None = None,
                 misfit=True, display_nnls=True, fail_gate=0.10, metric_cols=None, verbose=False,
                 se_method: str = "hc1") -> BranchFit:
    """Fit the branch engine for every HOG in `pairs` (polars/pandas DataFrame or dict of arrays).

    Required columns: hog, split, species_a, species_b, status and, per metric, value/se/sat columns
    (METRIC_COLS). `nu` = {metric: nu_m} tip correction (aa = 0). `kappa` = {metric: scalar or callable(d)}
    SE inflation (calibrated SE = kappa * Fisher SE). `hog_species` = {hog: S_h taxids} (default: species in pairs).
    """
    from .common import guard_split
    metric_cols = metric_cols or METRIC_COLS
    metrics = list(metrics)
    nu = {m: 0.0 for m in metrics} | (nu or {})
    nu["aa_ml"] = 0.0          # eta_aa = 0 (ANALYSIS_PLAN §1.4)
    kappa = kappa or {}
    hog = _col(pairs, "hog").astype(str)
    split = _col(pairs, "split").astype(str)
    guard_split(split)
    sa = _col(pairs, "species_a").astype(np.int64)
    sb = _col(pairs, "species_b").astype(np.int64)
    status = _col(pairs, "status").astype(str)
    vals, ses, oks = {}, {}, {}
    for m in metrics:
        vc, sc, fc = metric_cols[m]
        v = _col(pairs, vc).astype(float)
        s = _col(pairs, sc).astype(float)
        k = kappa.get(m, 1.0)
        s = s * (k(v) if callable(k) else k)
        sat = _col(pairs, fc).astype(bool) if fc is not None else np.zeros(len(v), bool)
        vals[m], ses[m] = v, s
        oks[m] = (status == "ok") & ~sat & np.isfinite(v) & np.isfinite(s)
    order = np.argsort(hog, kind="stable")
    hs = hog[order]
    bounds = np.flatnonzero(np.r_[True, hs[1:] != hs[:-1], True])
    K = len(metrics)
    aa_k = metrics.index("aa_ml") if "aa_ml" in metrics else None
    E = {k: [] for k in ("hog", "split", "metric", "edge_id", "is_terminal", "panel_child", "panel_parent",
                         "panel_edge", "n_panel_edges", "t_child", "t_parent", "tau", "b_raw", "b", "b_se",
                         "cov_aa", "b_nnls", "root_merged", "share", "_hog_code", "_child_idx", "_parent_idx",
                         "_param")}
    F = {k: [] for k in ("hog", "split", "metric", "n_pairs", "n_pairs_ok", "n_edges", "n_params", "s2_hm",
                         "rss_rel", "fail_frac", "misfit_ratio", "excluded", "reason", "var_sum", "cov_sum_aa",
                         "b_sum", "_hog_code")}
    P = {"row": [], "hidden_para": [], "resid_aa_z": []}
    hogfits = []
    for hc, (s0, s1) in enumerate(zip(bounds[:-1], bounds[1:])):
        rows = order[s0:s1]
        h = hs[s0]
        S = np.asarray(hog_species[h]) if hog_species is not None else np.unique(np.r_[sa[rows], sb[rows]])
        ht = prune(panel, S)
        pos = {int(t): k for k, t in enumerate(ht.species)}
        ia = np.array([pos[int(t)] for t in sa[rows]])
        ib = np.array([pos[int(t)] for t in sb[rows]])
        X = design_matrix(ht, ia, ib)
        d = np.stack([vals[m][rows] for m in metrics])
        se = np.stack([ses[m][rows] for m in metrics])
        ok = np.stack([oks[m][rows] for m in metrics])
        fit = wls_two_pass(X, d, se, ok, se_method=se_method)
        n = len(rows)
        fail = 1 - ok.mean(1)
        hf = HogFit(hog=h, split=split[rows[0]], ht=ht, ia=ia, ib=ib, X=X, metrics=metrics, d=d, se=se, ok=ok,
                    fit=fit, pair_index=rows)
        # misfit (AA, NJ topology) -- ANALYSIS_PLAN §2.5.6
        if misfit and aa_k is not None and not fit["under"][aa_k] and len(ht.species) >= 4:
            nsp = len(ht.species)
            fitted = X @ fit["b"][aa_k]
            D = np.zeros((nsp, nsp))
            dfill = np.where(ok[aa_k], d[aa_k], fitted)
            D[ia, ib] = dfill
            D[ib, ia] = dfill
            Mnj = nj_splits(D)
            Xnj = (Mnj[ia] ^ Mnj[ib]).astype(float)
            fnj = wls_two_pass(Xnj, d[aa_k:aa_k + 1], se[aa_k:aa_k + 1], ok[aa_k:aa_k + 1])
            if np.isfinite(fnj["rss_rel"][0]) and fnj["rss_rel"][0] > 0:
                hf.misfit_ratio = float(fit["rss_rel"][aa_k] / fnj["rss_rel"][0])
        # hidden paralogy (AA residual / 1.4826 MAD > 3)
        if aa_k is not None and not fit["under"][aa_k]:
            r = fit["r"][aa_k]
            okk = ok[aa_k]
            med = np.median(r[okk])
            mad = 1.4826 * np.median(np.abs(r[okk] - med))
            z = np.where(okk, r / mad if mad > 0 else 0.0, np.nan)
            hf.hidden_para = okk & (z > 3)
            P["row"].append(rows); P["hidden_para"].append(hf.hidden_para); P["resid_aa_z"].append(z)
        tc, tp, tau = edge_ages(panel, ht)
        eids = ht.edge_ids(panel)
        Ecount = len(ht.edge_child)
        for k, m in enumerate(metrics):
            reasons = []
            if fit["under"][k]:
                reasons.append("rank_deficient" if fit.get("rank_def", np.zeros(len(metrics), bool))[k] else "underdetermined")
            if fail[k] > fail_gate:
                reasons.append("fail_frac")
            if tau_misfit is not None and np.isfinite(hf.misfit_ratio) and hf.misfit_ratio > tau_misfit:
                reasons.append("misfit")
            hf.excluded[m] = ";".join(reasons)
            bp = fit["b"][k]
            Vp = fit["V"][k]
            if aa_k is not None and not fit["under"][aa_k] and not fit["under"][k]:
                if fit["se_method"] == "model":   # model-based has no residual cross-term: use HC1 scores
                    G1 = fit["A"][aa_k] * fit["r"][aa_k][None] * np.sqrt(fit["hc1"][aa_k])
                    G2 = fit["A"][k] * fit["r"][k][None] * np.sqrt(fit["hc1"][k])
                else:
                    G1, G2 = fit["G"][aa_k], fit["G"][k]
                Cp = G1 @ G2.T
            else:
                Cp = np.full_like(Vp, np.nan)
            bn = nnls_fit(X, d[k], fit["w"][k], ok[k]) if (display_nnls and not fit["under"][k]) else np.full(ht.n_params, np.nan)
            pe = ht.param_of_edge
            f = ht.share
            b_raw = f * bp[pe]
            b_se = np.abs(f) * np.sqrt(np.diag(Vp)[pe])
            cov = f * f * np.diag(Cp)[pe]
            b_corr = b_raw - np.where(ht.is_terminal, nu.get(m, 0.0), 0.0)
            # variance of the HOG total (all edges): 1'V1 on parameters (shares sum to 1 on merged params)
            ones = np.ones(ht.n_params)
            F["hog"].append(h); F["split"].append(hf.split); F["metric"].append(m); F["n_pairs"].append(n)
            F["n_pairs_ok"].append(int(fit["n_ok"][k])); F["n_edges"].append(Ecount); F["n_params"].append(ht.n_params)
            F["s2_hm"].append(fit["s2"][k]); F["rss_rel"].append(fit["rss_rel"][k]); F["fail_frac"].append(fail[k])
            F["misfit_ratio"].append(hf.misfit_ratio); F["excluded"].append(bool(reasons)); F["reason"].append(hf.excluded[m])
            F["var_sum"].append(float(ones @ Vp @ ones)); F["cov_sum_aa"].append(float(ones @ Cp @ ones))
            F["b_sum"].append(float(np.sum(b_corr))); F["_hog_code"].append(hc)
            for e in range(Ecount):
                c, pnode = ht.edge_child[e], ht.edge_parent[e]
                E["hog"].append(h); E["split"].append(hf.split); E["metric"].append(m); E["edge_id"].append(eids[e])
                E["is_terminal"].append(bool(ht.is_terminal[e])); E["panel_child"].append(panel.label[c])
                E["panel_parent"].append(panel.label[pnode])
                E["panel_edge"].append(panel.label[c] if ht.n_panel_edges[e] == 1 and f[e] == 1.0 else None)
                E["n_panel_edges"].append(int(ht.n_panel_edges[e])); E["t_child"].append(tc[e]); E["t_parent"].append(tp[e])
                E["tau"].append(tau[e]); E["b_raw"].append(b_raw[e]); E["b"].append(b_corr[e]); E["b_se"].append(b_se[e])
                E["cov_aa"].append(cov[e]); E["b_nnls"].append(f[e] * bn[pe[e]]); E["root_merged"].append(bool(f[e] != 1.0))
                E["share"].append(f[e]); E["_hog_code"].append(hc); E["_child_idx"].append(int(c))
                E["_parent_idx"].append(int(pnode)); E["_param"].append(int(pe[e]))
        hogfits.append(hf)
        if verbose and hc % 500 == 0:
            print(f"  branch engine: {hc} HOGs", flush=True)
    edges = {k: np.asarray(v) if k not in ("panel_edge",) else np.asarray(v, dtype=object) for k, v in E.items()}
    fits = {k: np.asarray(v) for k, v in F.items()}
    pf = {k: (np.concatenate(v) if v else np.zeros(0)) for k, v in P.items()}
    return BranchFit(aln=aln, metrics=metrics, nu=nu, hogs=hogfits, edges=edges, fits=fits, pair_flags=pf)


def refit_with(hf: HogFit, k: int, d_new: np.ndarray) -> np.ndarray:
    """Re-estimate parameters for new responses with the HOG's second-pass weights held fixed: b = A d.
    d_new (..., n). Returns parameters (..., p). Exact for any d_new (WLS is linear in d)."""
    return d_new @ hf.fit["A"][k].T


def edge_values_from_params(hf: HogFit, bp: np.ndarray, nu_m: float = 0.0) -> np.ndarray:
    """Map parameter vector(s) (..., p) to tip-corrected HOG-edge values (..., E)."""
    ht = hf.ht
    v = ht.share * bp[..., ht.param_of_edge]
    return v - np.where(ht.is_terminal, nu_m, 0.0)


def save_branch_tables(bf: BranchFit, panel: PanelTree, root=None, extra_fits=None):
    """Write data/interim/branches/{aln}.parquet (contract §1.3; contains ages -> interim only) and
    {aln}_fit.parquet. Also writes {aln}_pairflags.parquet (hidden-paralogy flags, row index into the pair table)."""
    from pathlib import Path
    import polars as pl
    from .common import ROOT
    out = Path(root or ROOT) / "data" / "interim" / "branches"
    out.mkdir(parents=True, exist_ok=True)
    bf.edge_frame(panel).write_parquet(out / f"{bf.aln}.parquet")
    bf.fit_frame().write_parquet(out / f"{bf.aln}_fit.parquet")
    pf = bf.pair_flags
    pl.DataFrame({"row": pf["row"].astype(np.int64), "hidden_para": pf["hidden_para"].astype(bool),
                  "resid_aa_z": pf["resid_aa_z"].astype(float)}).write_parquet(out / f"{bf.aln}_pairflags.parquet")
    return out
