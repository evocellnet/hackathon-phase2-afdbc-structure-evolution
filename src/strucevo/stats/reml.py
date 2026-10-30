"""Fast REML for the E4 mixed model  y = X beta + Z_h b_h + e,  b_h ~ N(0, sigma^2 Gamma),  e ~ N(0, sigma^2 I),
Z_h = [1, x] (random intercept + slope), small groups (<= ~10 rows) -- vectorised over groups by padding.

Same model and criterion as statsmodels MixedLM(reml=True) (profiled REML), but ~50-200x faster for many
small groups, which makes B_E4 = 1000 age-redraw bootstrap replicates affordable. Agreement with statsmodels
is checked at the point estimate in run_E4 and in tests/test_stats_engine.py (DEVIATION D6: used only inside
bootstrap replicates; the point estimate is statsmodels).
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize


def _pad(groups, *arrays):
    g = np.asarray(groups)
    u, inv = np.unique(g, return_inverse=True)
    order = np.argsort(inv, kind="stable")
    inv_s = inv[order]
    cnt = np.bincount(inv_s, minlength=len(u))
    m = int(cnt.max())
    start = np.r_[0, np.cumsum(cnt)[:-1]]
    pos = np.arange(len(inv_s)) - start[inv_s]
    out = []
    for a in arrays:
        a = np.asarray(a, float)
        shape = (len(u), m) + a.shape[1:]
        P = np.zeros(shape)
        P[inv_s, pos] = a[order]
        out.append(P)
    mask = np.zeros((len(u), m), bool)
    mask[inv_s, pos] = True
    return out, mask


class FastREML:
    def __init__(self, y, X, Z, groups):
        (self.Y, self.X, self.Z), self.mask = _pad(groups, y, X, Z)
        self.N = len(y)
        self.p = X.shape[1]
        self.q = Z.shape[1]
        self.m = self.Y.shape[1]
        self.I = np.eye(self.m)[None]

    def _L(self, th):
        L = np.zeros((self.q, self.q))
        L[np.tril_indices(self.q)] = th
        return L

    def crit(self, th):
        L = self._L(th)
        G = L @ L.T
        M = self.I + self.Z @ G @ np.swapaxes(self.Z, 1, 2)          # (H,m,m); padded rows: identity
        try:
            c = np.linalg.cholesky(M)
        except np.linalg.LinAlgError:
            return 1e30, None
        logdet = 2 * np.log(np.diagonal(c, axis1=1, axis2=2)).sum()
        RHS = np.concatenate([self.X, self.Y[:, :, None]], axis=2)    # (H,m,p+1)
        S = np.linalg.solve(M, RHS)
        XtMX = np.einsum("hmi,hmj->ij", self.X, S[:, :, :self.p])
        XtMy = np.einsum("hmi,hm->i", self.X, S[:, :, self.p])
        ytMy = np.einsum("hm,hm->", self.Y, S[:, :, self.p])
        try:
            beta = np.linalg.solve(XtMX, XtMy)
        except np.linalg.LinAlgError:
            return 1e30, None
        rss = ytMy - XtMy @ beta
        dof = self.N - self.p
        if rss <= 0:
            return 1e30, None
        s2 = rss / dof
        sign, ld2 = np.linalg.slogdet(XtMX)
        val = 0.5 * (dof * np.log(s2) + logdet + ld2)
        return val, (beta, s2, G, XtMX)

    def fit(self, starts=None):
        nq = self.q * (self.q + 1) // 2
        starts = starts or [np.r_[0.5, 0.0, 0.2][:nq] if nq == 3 else np.full(nq, 0.5), np.r_[1.5, 0.0, 0.5][:nq] if nq == 3 else np.full(nq, 1.0)]
        best = None
        for s in starts:
            r = minimize(lambda t: self.crit(t)[0], s, method="L-BFGS-B")
            if best is None or r.fun < best.fun:
                best = r
        val, (beta, s2, G, XtMX) = self.crit(best.x)
        cov_fe = s2 * np.linalg.inv(XtMX)
        return {"fe": beta, "s2": s2, "cov_re": s2 * G, "converged": bool(best.success), "crit": val,
                "se_fe": np.sqrt(np.diag(cov_fe)), "theta": best.x}


def fast_mixed(y, x, z, groups, slope=True, start=None):
    X = np.column_stack([np.ones_like(x), x, z])
    Z = np.column_stack([np.ones_like(x), x]) if slope else np.ones((len(x), 1))
    fr = FastREML(y, X, Z, groups)
    r = fr.fit(starts=[start] if start is not None else None)
    return {"mu": r["fe"][0], "beta": r["fe"][1], "gamma": r["fe"][2], "converged": r["converged"], "theta": r["theta"],
            "se_model": r["se_fe"]}
