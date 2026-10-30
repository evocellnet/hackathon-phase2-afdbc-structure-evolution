"""Pairwise curvature engine (ANALYSIS_PLAN §4.1/§4.2; DESIGN §11.1; R2-M3; R3-M6/M9).

Model  y = eta + a_h x^c + e  (eta fixed, a_h profiled per HOG, one global c).
  a_h(c) = sum_h w (y-eta) x^c / sum_h w x^{2c}
  Q(c)   = sum w (y-eta)^2 - sum_h [sum_h w (y-eta) x^c]^2 / sum_h w x^{2c}
  c_hat  = argmin Q: grid 0.20:0.01:3.00, then bounded scalar minimisation inside the grid bracket
           (scipy 'bounded' Brent for single fits; a vectorised golden-section search to 1e-9 for batched
           replicate fits -- identical minimiser, tested to 1e-6).
Two-pass weights: pass 1 w = 1/SE_y^2 ; s2 = max(0, mean r^2 - mean SE_y^2) ; pass 2 w = 1/(SE_y^2 + s2).
SIMEX: lambda in {0.5,1,1.5,2}, B_s replicates, bivariate N(0, lambda*Sigma_ij) with Var = kx^2 SE_x^2,
ky^2 SE_y^2, Cov = rho kx ky SE_x SE_y; x clipped at 1e-4 (clips counted); quadratic in lambda (incl. 0),
extrapolated to lambda = -1.  Weights are held at their pass-2 values inside SIMEX and bootstrap replicates
(DEVIATION/clarification D5).

All heavy evaluations run in torch float64 on the GPU when available (CPU otherwise).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import torch
from scipy.optimize import minimize_scalar

GRID = np.round(np.arange(0.20, 3.0 + 1e-9, 0.01), 2)
LAMBDAS = (0.5, 1.0, 1.5, 2.0)
GOLD = (math.sqrt(5) - 1) / 2


def _dev(device=None):
    if device is not None:
        return torch.device(device)
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


@dataclass
class CurveData:
    g: np.ndarray        # (n,) HOG code 0..H-1 (contiguous)
    x: np.ndarray
    y: np.ndarray
    sx: np.ndarray       # Fisher SEs of x (kappa applied inside SIMEX)
    sy: np.ndarray
    w0: np.ndarray       # eligibility (1/0) -- failed/saturated pairs carry 0
    eta: float
    H: int
    hog_ids: np.ndarray | None = None


class CurveEngine:
    """Vectorised evaluator of the profiled objective for many replicate data sets at once."""

    def __init__(self, data: CurveData, device=None, max_elems=6e7):
        self.dev = _dev(device)
        t = lambda a: torch.as_tensor(np.asarray(a, float), dtype=torch.float64, device=self.dev)
        self.d = data
        self.g = torch.as_tensor(np.asarray(data.g, np.int64), device=self.dev)
        self.x, self.y, self.sx, self.sy, self.w0 = t(data.x), t(data.y), t(data.sx), t(data.sy), t(data.w0)
        self.eta = float(data.eta)
        self.H = int(data.H)
        self.n = len(data.x)
        self.max_elems = max_elems

    # ---------------------------------------------------------------- core
    def _segsum(self, v):
        """v (..., n) -> (..., H) sums by HOG."""
        out = torch.zeros(v.shape[:-1] + (self.H,), dtype=v.dtype, device=self.dev)
        return out.index_add_(v.dim() - 1, self.g, v)

    def Q(self, c, x=None, y=None, w=None):
        """Profiled objective for replicate data. c: (R,G) or (G,) values; x, y, w: (R,n) or (n,) tensors.
        Returns (R,G)."""
        x = self.x if x is None else x
        y = self.y if y is None else y
        w = self.w if w is None else w
        if x.dim() == 1:
            x = x[None]
        if y.dim() == 1:
            y = y[None]
        if w.dim() == 1:
            w = w[None]
        R = max(x.shape[0], y.shape[0], w.shape[0])
        c = torch.as_tensor(c, dtype=torch.float64, device=self.dev)
        if c.dim() == 1:
            c = c[None].expand(R, -1)
        G = c.shape[1]
        lx = torch.log(x)
        ey = y - self.eta
        wy = w * ey                                                   # (R,n)
        base = (w * ey * ey).sum(-1)                                  # (R,)
        num = torch.zeros(R, G, self.H, dtype=torch.float64, device=self.dev)
        den = torch.zeros_like(num)
        # chunk over pairs so that R*G*chunk <= max_elems
        step = max(1, int(self.max_elems // max(1, R * G)))
        for s in range(0, self.n, step):
            e = min(self.n, s + step)
            lxs = lx[:, s:e] if lx.shape[0] > 1 else lx[:, s:e].expand(R, -1)
            xc = torch.exp(c[:, :, None] * lxs[:, None, :])           # (R,G,m)
            wys = (wy[:, s:e] if wy.shape[0] > 1 else wy[:, s:e].expand(R, -1))[:, None, :]
            ws = (w[:, s:e] if w.shape[0] > 1 else w[:, s:e].expand(R, -1))[:, None, :]
            gi = self.g[s:e]
            num.index_add_(2, gi, wys * xc)
            den.index_add_(2, gi, ws * xc * xc)
        prof = torch.where(den > 0, num * num / torch.where(den > 0, den, torch.ones_like(den)), torch.zeros_like(den))
        return base[:, None] - prof.sum(-1)

    def a_h(self, c, x=None, y=None, w=None):
        x = self.x if x is None else x
        y = self.y if y is None else y
        w = self.w if w is None else w
        xc = torch.exp(float(c) * torch.log(x))
        num = self._segsum(w * (y - self.eta) * xc)
        den = self._segsum(w * xc * xc)
        return torch.where(den > 0, num / torch.where(den > 0, den, torch.ones_like(den)), torch.full_like(den, float("nan")))

    def set_weights(self, w):
        self.w = torch.as_tensor(np.asarray(w, float), dtype=torch.float64, device=self.dev) if not torch.is_tensor(w) else w

    # ---------------------------------------------------------------- minimisers
    def fit_single(self, w=None, x=None, y=None, grid=GRID):
        """Grid + scipy bounded Brent. Returns (c_hat, at_boundary, Qmin)."""
        q = self.Q(torch.as_tensor(grid), x=x, y=y, w=w)[0].cpu().numpy()
        i = int(np.nanargmin(q))
        boundary = i == 0 or i == len(grid) - 1
        lo, hi = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
        f = lambda c: float(self.Q(torch.tensor([c]), x=x, y=y, w=w)[0, 0])
        r = minimize_scalar(f, bounds=(lo, hi), method="bounded", options={"xatol": 1e-9})
        return float(r.x), boundary, float(r.fun)


def golden_batch(engine: CurveEngine, lo, hi, x=None, y=None, w=None, tol=1e-9):
    """Clean vectorised golden-section search on brackets [lo, hi] (R,)."""
    dev = engine.dev
    a = torch.as_tensor(lo, dtype=torch.float64, device=dev).clone()
    b = torch.as_tensor(hi, dtype=torch.float64, device=dev).clone()
    Qc = lambda cc: engine.Q(cc[:, None], x=x, y=y, w=w)[:, 0]
    c = b - GOLD * (b - a)
    d = a + GOLD * (b - a)
    fc, fd = Qc(c), Qc(d)
    for _ in range(100):
        if float((b - a).max()) <= tol:
            break
        left = fc < fd                     # minimum in [a, d]
        # left: b<-d, d<-c, c<-new ; right: a<-c, c<-d, d<-new
        b = torch.where(left, d, b)
        a = torch.where(left, a, c)
        newpt = torch.where(left, b - GOLD * (b - a), a + GOLD * (b - a))
        fnew = Qc(newpt)
        c_new = torch.where(left, newpt, d)
        fc_new = torch.where(left, fnew, fd)
        d_new = torch.where(left, c, newpt)
        fd_new = torch.where(left, fc, fnew)
        c, d, fc, fd = c_new, d_new, fc_new, fd_new
    return (0.5 * (a + b)).cpu().numpy()


def fit_batch(engine: CurveEngine, x=None, y=None, w=None, grid=GRID, window=None, tol=1e-9):
    """Grid (optionally windowed, with automatic full-grid fallback) + vectorised golden section."""
    g_full = np.asarray(grid)
    g = g_full if window is None else g_full[(g_full >= window[0] - 1e-12) & (g_full <= window[1] + 1e-12)]
    q = engine.Q(torch.as_tensor(g), x=x, y=y, w=w)
    idx = torch.argmin(q, dim=1).cpu().numpy()
    R = q.shape[0]
    lo_c = g[np.maximum(idx - 1, 0)].copy()
    hi_c = g[np.minimum(idx + 1, len(g) - 1)].copy()
    boundary = (idx == 0) | (idx == len(g) - 1)
    if window is not None and boundary.any():
        sel = np.flatnonzero(boundary)

        def pick(a):
            if a is not None and torch.is_tensor(a) and a.dim() == 2 and a.shape[0] == R:
                return a[torch.as_tensor(sel, device=a.device)]
            return a
        qf = engine.Q(torch.as_tensor(g_full), x=pick(x), y=pick(y), w=pick(w))
        jf = torch.argmin(qf, dim=1).cpu().numpy()
        lo_c[sel] = g_full[np.maximum(jf - 1, 0)]
        hi_c[sel] = g_full[np.minimum(jf + 1, len(g_full) - 1)]
        boundary[sel] = (jf == 0) | (jf == len(g_full) - 1)
    return golden_batch(engine, lo_c, hi_c, x=x, y=y, w=w, tol=tol), boundary


# ------------------------------------------------------------------------------------------ fits
def two_pass(engine: CurveEngine, wmult=None, fixed_s2=None):
    """Two-pass weighted fit. wmult: bootstrap multiplicity weights (n,) or None.
    Returns dict(c, boundary, s2, w (tensor), a_h (np), Qmin)."""
    sy2 = engine.sy ** 2
    wm = engine.w0 if wmult is None else engine.w0 * torch.as_tensor(np.asarray(wmult, float), dtype=torch.float64, device=engine.dev)
    elig = engine.w0 > 0
    if fixed_s2 is None:
        w1 = torch.where(elig, wm / torch.where(elig, sy2, torch.ones_like(sy2)), torch.zeros_like(sy2))
        c1, _, _ = engine.fit_single(w=w1)
        ah = engine.a_h(c1, w=w1)
        fitted = engine.eta + ah[engine.g] * torch.exp(c1 * torch.log(engine.x))
        r = (engine.y - fitted)[elig]
        s2 = max(0.0, float((r * r).mean() - sy2[elig].mean()))
    else:
        s2 = float(fixed_s2)
    w2 = torch.where(elig, wm / (sy2 + s2), torch.zeros_like(sy2))
    c, boundary, qmin = engine.fit_single(w=w2)
    return {"c": c, "boundary": boundary, "s2": s2, "w": w2, "a_h": engine.a_h(c, w=w2).cpu().numpy(), "Qmin": qmin}


@dataclass
class TreeX:
    """Error structure of the tree-smoothed regressor x_hat = X_h b_hat_AA (DEVIATION D29).
    L: (H, p, p) square roots of V_h = Cov(b_hat_AA) (padded with zeros); Xrow: (n, p) path-incidence rows of the
    curve pairs (padded); sd: (n,) SD of x_hat; rho_scale: (n,) Corr(x_hat error, own AA error) so that the y noise
    gets Corr(e_x_hat, e_y) = rho * rho_scale (own-pair term of the hat matrix)."""
    L: object
    Xrow: object
    sd: object
    rho_scale: object


def simex(engine: CurveEngine, w, c0, kx=1.0, ky=1.0, rho=0.0, lambdas=LAMBDAS, B_s=50, rng=None,
          window=0.3, batch=None, xtree: TreeX | None = None, clip_mode="clip"):
    """SIMEX for c (quadratic extrapolant to lambda=-1). w: pass-2 weights (tensor, incl. bootstrap mult).
    kx, ky, rho: scalars or per-pair arrays (n,) (per-pair Sigma_ij, ANALYSIS_PLAN §4.1).
    xtree: if given, x is the tree-smoothed AA distance and its remeasurement noise is drawn with its exact
    within-HOG covariance X_h V_h X_h' (b* = b_hat + sqrt(lambda) L_h zeta); sx/kx are then ignored.
    Returns dict(c_simex, c_lambda, n_clip, boundary_any)."""
    rng = np.random.default_rng(rng)
    dev = engine.dev
    n = engine.n
    _t = lambda v: (torch.as_tensor(np.asarray(v, float), dtype=torch.float64, device=dev)
                    if np.ndim(v) else float(v))      # scalars stay scalars; per-pair arrays (n,) -> tensors
    kx, ky, rho = _t(kx), _t(ky), _t(rho)
    sx = engine.sx * kx
    sy = engine.sy * ky
    rho_c = (torch.sqrt(torch.clamp(1 - rho * rho, min=0.0)) if torch.is_tensor(rho)
             else math.sqrt(max(0.0, 1 - rho * rho)))
    win = None if window is None else (max(GRID[0], c0 - window), min(GRID[-1], c0 + window))
    c_lam, n_clip, bnd = [], 0, False
    if xtree is not None:
        p = xtree.L.shape[-1]
        batch = batch or max(1, min(B_s, int(6e7 // max(1, n * p))))
        rho_t = rho * xtree.rho_scale
        rho_tc = torch.sqrt(torch.clamp(1 - rho_t * rho_t, min=0.0))
    batch = batch or B_s
    for lam in lambdas:
        vals = []
        for s in range(0, B_s, batch):
            R = min(batch, B_s - s)
            z2 = torch.as_tensor(rng.standard_normal((R, n)), dtype=torch.float64, device=dev)
            if xtree is None:
                z1 = torch.as_tensor(rng.standard_normal((R, n)), dtype=torch.float64, device=dev)
                ex = math.sqrt(lam) * sx * z1
                ey = math.sqrt(lam) * sy * (rho * z1 + rho_c * z2)
            else:
                zeta = torch.as_tensor(rng.standard_normal((R,) + tuple(xtree.L.shape[:2])), dtype=torch.float64, device=dev)
                u = torch.einsum("hkl,rhl->rhk", xtree.L, zeta)                   # (R,H,p) edge-parameter errors
                z1 = (xtree.Xrow[None] * u[:, engine.g, :]).sum(-1) / xtree.sd   # standardised x_hat errors
                ex = math.sqrt(lam) * xtree.sd * z1
                ey = math.sqrt(lam) * sy * (rho_t * z1 + rho_tc * z2)
            xs = engine.x + ex
            clip = xs < 1e-4
            n_clip += int((clip & (engine.w0 > 0)).sum())
            xs = torch.clamp(xs, min=1e-4)
            ys = engine.y + ey
            wr = w if clip_mode == "clip" else torch.where(clip, torch.zeros_like(xs), w if w.dim() == 2 else w[None].expand(R, -1))
            c, b = fit_batch(engine, x=xs, y=ys, w=wr, window=win)
            bnd |= bool(b.any())
            vals.append(c)
        c_lam.append(float(np.mean(np.concatenate(vals))))
    lam_all = np.r_[0.0, lambdas]
    cs = np.r_[c0, c_lam]
    coef = np.polyfit(lam_all, cs, 2)
    return {"c_simex": float(np.polyval(coef, -1.0)), "c_lambda": dict(zip(map(float, lam_all), map(float, cs))),
            "n_clip": n_clip, "boundary_any": bnd, "quad_coef": coef.tolist()}


def weighted_median(v, w):
    v = np.asarray(v, float); w = np.asarray(w, float)
    m = np.isfinite(v) & np.isfinite(w) & (w > 0)
    v, w = v[m], w[m]
    if v.size == 0:
        return float("nan")
    o = np.argsort(v)
    cw = np.cumsum(w[o])
    return float(v[o][np.searchsorted(cw, 0.5 * cw[-1])])


def fit_eta_free(engine: CurveEngine, w, eta_bounds=(-0.5, 0.5)):
    """Sensitivity 'eta free': minimise min_c Q(c; eta) over eta (bounded scalar)."""
    eta0 = engine.eta

    def prof(eta):
        engine.eta = float(eta)
        c, _, q = engine.fit_single(w=w)
        return q
    r = minimize_scalar(prof, bounds=(eta0 + eta_bounds[0], eta0 + eta_bounds[1]), method="bounded",
                        options={"xatol": 1e-6})
    engine.eta = float(r.x)
    c, b, q = engine.fit_single(w=w)
    engine.eta = eta0
    return {"eta": float(r.x), "c": c, "boundary": b}


# ------------------------------------------------------------------------------------------ pure numpy reference
def reference_fit(g, x, y, w, eta, grid=GRID):
    """Independent CPU re-implementation (for the 1e-6 recomputation check)."""
    g = np.asarray(g); H = g.max() + 1

    def Q(c):
        xc = x ** c
        num = np.bincount(g, w * (y - eta) * xc, H)
        den = np.bincount(g, w * xc * xc, H)
        ok = den > 0
        return np.sum(w * (y - eta) ** 2) - np.sum(num[ok] ** 2 / den[ok])
    q = np.array([Q(c) for c in grid])
    i = int(np.argmin(q))
    r = minimize_scalar(Q, bounds=(grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]), method="bounded",
                        options={"xatol": 1e-10})
    return float(r.x)


# ------------------------------------------------------------------------------------------ D29: SIMEX bias correction
def simex_pb_correction(engine: CurveEngine, fit2: dict, c_simex: float, kx=1.0, ky=1.0, rho=0.0, K=25, B_s=20,
                        rng=None, lambdas=LAMBDAS, window=0.3, x_design=None):
    """DEVIATION D29 (agents/bias_fix_report.md): parametric-bootstrap correction of the residual SIMEX bias.

    The quadratic SIMEX extrapolant under-corrects when the naive attenuation is large (realistic AA SEs). Pseudo
    worlds take the fitted model as the truth: D* = observed x (design), y* mean = eta + a_h(c_S) x^c_S with
    c_S = the SIMEX estimate; observation noise = the calibrated per-pair Sigma_ij (same kx, ky, rho as SIMEX) plus
    N(0, s2) on y. Each pseudo data set is re-estimated with the identical estimator (two-pass + SIMEX, B_s).
    bias = mean(c*_S) - c_S ; c_corrected = c_S - bias. Returns dict(c_corrected, bias, mc_se, reps).
    The correction is computed once on the full data and held fixed in bootstrap replicates (as the pass-2
    weights, D5); its Monte-Carlo SE is added in quadrature to the bootstrap SE by the caller."""
    rng = np.random.default_rng(rng)
    dev = engine.dev
    d = engine.d
    _t = lambda v: (torch.as_tensor(np.asarray(v, float), dtype=torch.float64, device=dev) if np.ndim(v) else float(v))
    kx_, ky_, rho_ = _t(kx), _t(ky), _t(rho)
    sx = engine.sx * kx_
    sy = engine.sy * ky_
    rho_c = torch.sqrt(torch.clamp(1 - rho_ * rho_, min=0.0)) if torch.is_tensor(rho_) else math.sqrt(max(0.0, 1 - rho_ ** 2))
    cS = float(min(max(c_simex, GRID[0]), GRID[-1]))
    ah = engine.a_h(cS, w=fit2["w"])
    xd = engine.x if x_design is None else torch.as_tensor(np.asarray(x_design, float), dtype=torch.float64, device=dev)
    mu = engine.eta + ah[engine.g] * torch.exp(cS * torch.log(xd))
    s2 = float(fit2["s2"])
    elig = engine.w0 > 0
    reps = []
    for k in range(K):
        z1 = torch.as_tensor(rng.standard_normal(engine.n), dtype=torch.float64, device=dev)
        z2 = torch.as_tensor(rng.standard_normal(engine.n), dtype=torch.float64, device=dev)
        z3 = torch.as_tensor(rng.standard_normal(engine.n), dtype=torch.float64, device=dev)
        xs = torch.clamp(xd + sx * z1, min=1e-4)
        ys = mu + sy * (rho_ * z1 + rho_c * z2) + math.sqrt(s2) * z3
        xs = torch.where(elig, xs, engine.x); ys = torch.where(elig, ys, engine.y)
        dd = CurveData(g=d.g, x=xs.cpu().numpy(), y=ys.cpu().numpy(), sx=d.sx, sy=d.sy, w0=d.w0, eta=engine.eta, H=d.H)
        e2 = CurveEngine(dd, device=dev)
        f2 = two_pass(e2)
        sm = simex(e2, f2["w"], f2["c"], kx=kx, ky=ky, rho=rho, lambdas=lambdas, B_s=B_s, rng=rng, window=window)
        reps.append(sm["c_simex"])
    reps = np.asarray(reps)
    bias = float(reps.mean() - cS)
    return {"c_corrected": float(c_simex - bias), "bias": bias, "mc_se": float(reps.std(ddof=1) / math.sqrt(K)),
            "reps": reps, "K": K, "B_s": B_s}
