"""Curve models for divergence y versus x (e.g. 3Di / TM distance vs AA distance or time).

Models
    lin0  y = b x                    (linear through origin)            k = 1
    lin   y = a + b x                (linear with intercept)            k = 2
    sat   y = a (1 - exp(-b x))      (saturating exponential)           k = 2
    power y = a x^c                  (power law, x > 0)                 k = 2
All fitted by (weighted) least squares. Nonlinear models use variable projection: for each value of the
nonlinear parameter on a wide log grid the linear amplitude has a closed form; the best grid point is then
refined with scipy.optimize.least_squares -> robust to starting values.

Information criteria assume i.i.d. Gaussian errors with unknown variance (counted as a parameter):
    AIC = n ln(SSE/n) + 2(k+1),  AICc = AIC + 2(k+1)(k+2)/(n-k-2),  BIC = n ln(SSE/n) + (k+1) ln n.
Out-of-range error: fit on x < quantile(x, q), predict x >= quantile -> RMSE, MAE, mean signed error.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
from scipy.optimize import least_squares

MODEL_NAMES = ("lin0", "lin", "sat", "power")


@dataclass
class FitResult:
    name: str
    params: dict
    sse: float
    n: int
    k: int
    aic: float
    aicc: float
    bic: float
    predict: Callable = field(repr=False)
    success: bool = True

    def as_dict(self):
        return dict(name=self.name, **self.params, sse=self.sse, n=self.n, k=self.k, aic=self.aic,
                    aicc=self.aicc, bic=self.bic, success=self.success)


def _f(name, p, x):
    x = np.asarray(x, float)
    if name == "lin0":
        return p["b"] * x
    if name == "lin":
        return p["a"] + p["b"] * x
    if name == "sat":
        return p["a"] * (-np.expm1(-p["b"] * x))
    if name == "power":
        with np.errstate(invalid="ignore", divide="ignore"):
            return p["a"] * np.power(np.clip(x, 0, None), p["c"])
    raise KeyError(name)


def _ic(sse, n, k):
    sse = max(sse, 1e-300)
    kk = k + 1
    aic = n * np.log(sse / n) + 2 * kk
    aicc = aic + (2 * kk * (kk + 1) / (n - kk - 1) if n - kk - 1 > 0 else np.inf)
    bic = n * np.log(sse / n) + kk * np.log(n)
    return aic, aicc, bic


def _amp(g, y, w):
    """Closed-form LS amplitude a for y ~ a*g."""
    den = (w * g * g).sum()
    return (w * g * y).sum() / den if den > 0 else 0.0


def fit_model(name: str, x, y, weights=None) -> FitResult:
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    w = np.ones_like(x) if weights is None else np.asarray(weights, float)
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
    if name == "power":
        ok &= x > 0
    x, y, w = x[ok], y[ok], w[ok]
    n = len(x)
    sw = np.sqrt(w)
    success = True
    if name == "lin0":
        p = dict(b=_amp(x, y, w))
        k = 1
    elif name == "lin":
        X = np.c_[np.ones(n), x] * sw[:, None]
        coef, *_ = np.linalg.lstsq(X, y * sw, rcond=None)
        p = dict(a=coef[0], b=coef[1])
        k = 2
    elif name == "sat":
        k = 2
        xs = np.median(x[x > 0]) if np.any(x > 0) else 1.0
        bgrid = np.logspace(-4, 3, 300) / xs
        best = None
        for b in bgrid:
            g = -np.expm1(-b * x)
            a = _amp(g, y, w)
            sse = (w * (y - a * g) ** 2).sum()
            if best is None or sse < best[0]:
                best = (sse, a, b)
        _, a0, b0 = best
        res = least_squares(lambda t: sw * (y - t[0] * (-np.expm1(-np.exp(t[1]) * x))),
                            x0=[a0, np.log(b0)], method="lm" if n > 2 else "trf")
        success = bool(res.success)
        p = dict(a=res.x[0], b=float(np.exp(res.x[1])))
        if (w * (y - _f("sat", p, x)) ** 2).sum() > best[0]:
            p = dict(a=a0, b=b0)
    elif name == "power":
        k = 2
        best = None
        for c in np.linspace(0.02, 4.0, 400):
            g = x ** c
            a = _amp(g, y, w)
            sse = (w * (y - a * g) ** 2).sum()
            if best is None or sse < best[0]:
                best = (sse, a, c)
        _, a0, c0 = best
        res = least_squares(lambda t: sw * (y - t[0] * x ** t[1]), x0=[a0, c0],
                            method="lm" if n > 2 else "trf")
        success = bool(res.success)
        p = dict(a=res.x[0], c=res.x[1])
        if not np.isfinite(p["c"]) or (w * (y - _f("power", p, x)) ** 2).sum() > best[0]:
            p = dict(a=a0, c=c0)
    else:
        raise KeyError(name)
    p = {kk: float(v) for kk, v in p.items()}
    sse = float((w * (y - _f(name, p, x)) ** 2).sum())
    aic, aicc, bic = _ic(sse, n, k)
    return FitResult(name, p, sse, n, k, float(aic), float(aicc), float(bic),
                     predict=lambda xx, _n=name, _p=dict(p): _f(_n, _p, xx), success=success)


def fit_models(x, y, models=MODEL_NAMES, weights=None) -> dict:
    """Fit several models to the same data -> {name: FitResult}.

    ICs are only comparable on identical data: if 'power' is among the models, points with x <= 0 are dropped
    for ALL models (power is undefined there).
    """
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if "power" in models:
        keep = x > 0
        x, y = x[keep], y[keep]
        weights = None if weights is None else np.asarray(weights, float)[keep]
    return {m: fit_model(m, x, y, weights) for m in models}


def compare(fits: dict, criterion: str = "aic") -> list:
    """Table (list of dicts) sorted by criterion with delta and Akaike/Schwarz weights."""
    vals = np.array([getattr(f, criterion) for f in fits.values()])
    delta = vals - vals.min()
    wts = np.exp(-0.5 * delta)
    wts /= wts.sum()
    rows = []
    for (name, f), dl, wt in zip(fits.items(), delta, wts):
        r = f.as_dict()
        r[f"delta_{criterion}"] = float(dl)
        r[f"weight_{criterion}"] = float(wt)
        rows.append(r)
    return sorted(rows, key=lambda r: r[f"delta_{criterion}"])


def best_model(x, y, models=MODEL_NAMES, criterion: str = "aic", weights=None) -> str:
    fits = fit_models(x, y, models, weights)
    return min(fits, key=lambda m: getattr(fits[m], criterion))


def out_of_range_error(x, y, q: float = 0.5, models=MODEL_NAMES, weights=None) -> dict:
    """Fit on x < quantile(x, q), predict the rest. Returns {model: dict(rmse, mae, bias, n_train, n_test, cut)}."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if "power" in models:
        keep = x > 0
        x, y = x[keep], y[keep]
        weights = None if weights is None else np.asarray(weights, float)[keep]
    cut = float(np.quantile(x[np.isfinite(x)], q))
    tr, te = x < cut, x >= cut
    w = None if weights is None else np.asarray(weights, float)
    out = {}
    for m in models:
        f = fit_model(m, x[tr], y[tr], None if w is None else w[tr])
        r = y[te] - f.predict(x[te])
        out[m] = dict(rmse=float(np.sqrt(np.mean(r ** 2))), mae=float(np.mean(np.abs(r))),
                      bias=float(np.mean(-r)), n_train=int(tr.sum()), n_test=int(te.sum()), cut=cut,
                      params=f.params)
    return out


def divergence_x(f_a: FitResult, f_b: FitResult, rel: float = 0.10, x_max: float = 10.0, n: int = 10001):
    """Smallest x > 0 where |f_a(x) - f_b(x)| > rel * |f_b(x)| (e.g. linear vs saturating: d*); NaN if none."""
    xs = np.linspace(x_max / n, x_max, n)
    ya, yb = f_a.predict(xs), f_b.predict(xs)
    with np.errstate(invalid="ignore", divide="ignore"):
        hit = np.abs(ya - yb) > rel * np.abs(yb)
    return float(xs[np.argmax(hit)]) if hit.any() else float("nan")
