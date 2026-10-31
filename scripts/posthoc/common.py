"""Shared helpers for scripts/posthoc (copied from pairwise_checks.py): E1 window extraction and WLS curvature fits
with fixed eta / HOG intercepts / per-protein additive effects profiled out (Frisch-Waugh)."""
import numpy as np, scipy.sparse as sp
from scipy.optimize import minimize_scalar
from strucevo.stats.endpoints_curve import prepare_curve_data

def get_window(tab, excl, lo, hi):
    data, info, idx, _ = prepare_curve_data(tab, "E1", excluded_hogs=excl, x_range=(lo, hi))
    m = data.w0 > 0
    ii = idx[m]
    ga = tab["gene_a"].to_numpy()[ii].astype(str); gb = tab["gene_b"].to_numpy()[ii].astype(str)
    hogs = data.g[m]
    ug, g = np.unique(hogs, return_inverse=True)
    # protein ids unique within HOG
    ka = np.array([f"{a}|{b}" for a, b in zip(g, ga)]); kb = np.array([f"{a}|{b}" for a, b in zip(g, gb)])
    up, pinv = np.unique(np.r_[ka, kb], return_inverse=True)
    pa, pb = pinv[:len(ka)], pinv[len(ka):]
    phog = np.zeros(len(up), int); phog[pa] = g; phog[pb] = g
    return dict(x=data.x[m], y=data.y[m], sx=data.sx[m], sy=data.sy[m], g=g, H=len(ug), pa=pa, pb=pb, P=len(up),
                phog=phog, rows=ii, prot=up)

class Fitter:
    """Weighted least squares with nuisance blocks profiled out (Frisch-Waugh)."""
    def __init__(self, D, w, mode):
        self.D, self.w, self.mode = D, w, mode
        n = len(w)
        if mode == "hog":
            M = sp.csr_matrix((np.ones(n), (np.arange(n), D["g"])), shape=(n, D["H"]))
        elif mode == "prot":
            M = sp.csr_matrix((np.ones(2 * n), (np.r_[np.arange(n), np.arange(n)], np.r_[D["pa"], D["pb"]])), shape=(n, D["P"]))
        else:
            M = None
        self.M = M
        if M is not None:
            A = (M.T @ sp.diags(w) @ M).tocsr()
            # block-diagonal pseudo-inverse per HOG
            if mode == "hog":
                self.Ainv = sp.diags(1 / np.maximum(A.diagonal(), 1e-300))
            else:
                blocks = []
                order = np.argsort(D["phog"], kind="stable")
                assert np.all(order == np.arange(D["P"])) or True
                # proteins are sorted by "g|gene" string -> sort by g numerically to build blocks
                perm = np.argsort(D["phog"], kind="stable")
                self.perm = perm
                Ap = A[perm][:, perm].tocsr()
                bounds = np.r_[0, np.cumsum(np.bincount(D["phog"][perm], minlength=D["H"]))]
                for s, e in zip(bounds[:-1], bounds[1:]):
                    blocks.append(np.linalg.pinv(Ap[s:e, s:e].toarray(), rcond=1e-10, hermitian=True))
                Bi = sp.block_diag(blocks, format="csr")
                Pm = sp.csr_matrix((np.ones(D["P"]), (perm, np.arange(D["P"]))), shape=(D["P"], D["P"]))
                self.Ainv = (Pm @ Bi @ Pm.T).tocsr()
            self.MtW = (M.T @ sp.diags(w)).tocsr()

    def resid(self, v):
        if self.M is None:
            return v
        return v - self.M @ (self.Ainv @ (self.MtW @ v))

    def Q(self, c, y, eta=None):
        D, w = self.D, self.w
        u = D["xx"] ** c
        if self.mode == "fixed":
            yr = y - eta; ur = u
        else:
            yr = self.resid(y); ur = self.resid(u)
        g, H = D["g"], D["H"]
        suy = np.bincount(g, w * ur * yr, H); suu = np.bincount(g, w * ur * ur, H)
        a = suy / np.maximum(suu, 1e-300)
        r = yr - a[g] * ur
        return float(np.sum(w * r * r)), a, r

    def fit_c(self, y, eta=None):
        grid = np.arange(0.3, 1.81, 0.05)
        q = [self.Q(c, y, eta)[0] for c in grid]; i = int(np.argmin(q))
        o = minimize_scalar(lambda c: self.Q(c, y, eta)[0], bounds=(grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]),
                            method="bounded", options={"xatol": 1e-4})
        return float(o.x), float(o.fun), bool(i in (0, len(grid) - 1))

    def offsets(self, c, y, eta=None):
        """nuisance effects at curvature c."""
        u = self.D["xx"] ** c
        _, a, _ = self.Q(c, y, eta)
        resid_a = y - a[self.D["g"]] * u
        return self.Ainv @ (self.MtW @ resid_a)

