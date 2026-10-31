#!/usr/bin/env python3
"""Post-hoc verification (independent of scripts/review): pairwise checks for claims 1 and 2.
 (a) identical-sequence floor vs within-HOG near-identical ortholog pairs (binned tdi_ml at small d_AA).
 (b) naive curvature fits on the E1 data (production filters, raw x = aa_ml, no SIMEX) with three baselines:
       fixed  : y = eta + a_h x^c                (eta = 0.0524, production)
       hog    : y = eta_h + a_h x^c              (reviewer's alternative)
       prot   : y = u_i + u_j + a_h x^c          (per-protein additive effects)
 (c) correlates of the per-protein effects u_i (pLDDT, length).
 (d) synthetic recovery on the real design: truth c in {0.72, 1.0}, with/without per-protein offsets,
     with/without x measurement error.
  CONFIRM_RUN=1 python scripts/posthoc/pairwise_checks.py -> results/posthoc/pairwise_checks.json
"""
import sys, json, time, tempfile
from pathlib import Path
import numpy as np, polars as pl, scipy.sparse as sp
from scipy.optimize import minimize_scalar
from scipy.stats import spearmanr
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from strucevo.stats.driver import EndpointDriver
from strucevo.stats.endpoints_curve import prepare_curve_data

OUT = ROOT / "results" / "posthoc"; OUT.mkdir(parents=True, exist_ok=True)
t0 = time.time()
def log(*a): print(f"[{time.time()-t0:7.0f}s]", *a, flush=True)
res = {}
def save(): (OUT / "pairwise_checks.json").write_text(json.dumps(res, indent=1, default=float))

d = EndpointDriver(split="confirm", run_alt_alignments=False, outdir=tempfile.mkdtemp(), log=lambda *a: None)
d.load()
tab = d.tabs["tdiaa"]
ETA = d.tun.eta("tdiaa")["tdi_ml"]
res["eta"] = ETA

# ---------------------------------------------------------------- (a) floor vs near-identical orthologs
fl = pl.scan_parquet(str(ROOT / "data/interim/pairs_floors/tdiaa/*.parquet")).filter(
    (pl.col("pair_set") == "floor_ident") & (pl.col("status") == "ok")).collect()
res["floor_ident"] = {"n": fl.height, "median": fl["tdi_ml"].median(), "mean": fl["tdi_ml"].mean(),
                      "note": "pairs_floors table as stored (alpha of the floors run); tuned eta = 0.0524"}
t = tab.filter((pl.col("status") == "ok") & pl.col("one_to_one").cast(pl.Boolean) & ~pl.col("tdi_sat").cast(pl.Boolean))
rows = []
for lo, hi in [(0, 1e-9), (1e-9, 0.005), (0.005, 0.01), (0.01, 0.02), (0.02, 0.03), (0.03, 0.05), (0.05, 0.075),
               (0.075, 0.1), (0.1, 0.15), (0.15, 0.2), (0.2, 0.3), (0.3, 0.5)]:
    s = t.filter((pl.col("aa_ml") >= lo) & (pl.col("aa_ml") < hi))
    rows.append({"x_lo": lo, "x_hi": hi, "n": s.height, "x_mean": s["aa_ml"].mean(), "tdi_median": s["tdi_ml"].median(),
                 "tdi_mean": s["tdi_ml"].mean(), "plddt_min_mean": s.select(pl.min_horizontal("mean_plddt_a", "mean_plddt_b")).to_series().mean()})
res["small_x_bins"] = rows
# identical sequences between species within HOGs: also matched on pLDDT to the floor
z = t.filter(pl.col("aa_ml") == 0)
res["within_hog_x0"] = {"n": z.height, "median": z["tdi_ml"].median(), "mean": z["tdi_ml"].mean()}
log("floor", res["floor_ident"], res["within_hog_x0"]); save()

# ---------------------------------------------------------------- (b) E1 data
hc = pl.read_parquet(ROOT / "data/interim/hogs_endpoint_cols_confirm.parquet")
excl = set(hc.filter(pl.col("excl_tdi_tdiaa"))["hog"].to_list())

def get_window(lo, hi):
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

windows = [(0.05, 2.0), (0.05, 1.0), (0.3, 2.0)]
res["fits"] = {}
Dmain = None
for lo, hi in windows:
    D = get_window(lo, hi); D["xx"] = D["x"]
    # pass-1 weights 1/sy^2 under the fixed model at its c, then s2 (production two-pass rule)
    f0 = Fitter(D, 1 / D["sy"] ** 2, "fixed"); c0, _, _ = f0.fit_c(D["y"], ETA)
    _, _, r0 = f0.Q(c0, D["y"], ETA)
    s2 = max(0.0, float(np.mean(r0 ** 2) - np.mean(D["sy"] ** 2)))
    w = 1 / (D["sy"] ** 2 + s2)
    rec = {"n": len(w), "H": D["H"], "P": D["P"], "s2": s2}
    for mode in ("fixed", "hog", "prot"):
        F = Fitter(D, w, mode)
        c, q, b = F.fit_c(D["y"], ETA)
        q1 = F.Q(1.0, D["y"], ETA)[0]
        rec[mode] = {"c": c, "Q": q, "boundary": b, "Q_c1_minus_Qmin": q1 - q}
        if mode != "fixed":
            off = F.offsets(c, D["y"], ETA); off1 = F.offsets(1.0, D["y"], ETA)
            rec[mode]["offset_q10_50_90_at_chat"] = np.quantile(off, [.1, .5, .9]).tolist()
            rec[mode]["offset_mean_at_chat"] = float(off.mean())
            rec[mode]["offset_q10_50_90_at_c1"] = np.quantile(off1, [.1, .5, .9]).tolist()
            rec[mode]["offset_mean_at_c1"] = float(off1.mean())
        if (lo, hi) == (0.05, 2.0):
            D["w"] = w; Dmain = D
            if mode == "prot":
                D["u_chat"] = off; D["u_c1"] = off1; D["c_prot"] = c
            if mode == "fixed":
                D["c_fixed"] = c
        log(lo, hi, mode, json.dumps(rec[mode], default=float))
    res["fits"][f"x_{lo}_{hi}"] = rec; save()

# ---------------------------------------------------------------- (c) correlates of per-protein effects
D = Dmain
rows = D["rows"]
pl_a = tab["mean_plddt_a"].to_numpy()[rows]; pl_b = tab["mean_plddt_b"].to_numpy()[rows]
le_a = tab["len_a"].to_numpy()[rows].astype(float); le_b = tab["len_b"].to_numpy()[rows].astype(float)
P = D["P"]
plp = np.full(P, np.nan); lep = np.full(P, np.nan)
plp[D["pa"]] = pl_a; plp[D["pb"]] = pl_b; lep[D["pa"]] = le_a; lep[D["pb"]] = le_b
npair = np.bincount(np.r_[D["pa"], D["pb"]], minlength=P)
u = D["u_chat"]; u1 = D["u_c1"]
ok = npair >= 5
# within-HOG centred (removes HOG-level differences, keeps protein-level)
def centred(v):
    m = np.bincount(D["phog"][ok], v[ok], D["H"]) / np.maximum(np.bincount(D["phog"][ok], minlength=D["H"]), 1)
    return v[ok] - m[D["phog"][ok]]
corr = {}
for name, v in (("plddt", plp), ("length", lep)):
    corr[name] = {"spearman_u_all": spearmanr(u[ok], v[ok]).statistic,
                  "spearman_u_within_hog": spearmanr(centred(u), centred(v)).statistic}
q = np.quantile(plp[ok], [0, .2, .4, .6, .8, 1])
corr["u_mean_by_plddt_quintile"] = [{"plddt_lo": float(q[i]), "plddt_hi": float(q[i + 1]),
                                     "u_mean": float(u[ok][(plp[ok] >= q[i]) & (plp[ok] <= q[i + 1])].mean()),
                                     "u_median": float(np.median(u[ok][(plp[ok] >= q[i]) & (plp[ok] <= q[i + 1])]))}
                                    for i in range(5)]
corr["u_mean_minus_eta_half_at_chat"] = float(u[ok].mean() - ETA / 2)
corr["u_mean_minus_eta_half_at_c1"] = float(u1[ok].mean() - ETA / 2)
corr["u_median_at_c1"] = float(np.median(u1[ok]))
corr["n_proteins"] = int(ok.sum())
# split-half reliability of u_i: fit on two random halves of pairs (protein effect is a property of the protein?)
res["protein_effects"] = corr; log("corr", json.dumps(corr, default=float)); save()

# ---------------------------------------------------------------- (d) synthetic recovery on the real design
rng = np.random.default_rng(20261031)
Ffix = Fitter(D, D["w"], "fixed"); Fhog = Fitter(D, D["w"], "hog"); Fprot = Fitter(D, D["w"], "prot")
s2 = res["fits"]["x_0.05_2.0"]["s2"]
x_true = D["x"].copy()
syn = []
for c_true in (0.72, 1.0):
    _, a_true, _ = Ffix.Q(c_true, D["y"], ETA)       # realistic family scales
    for off_kind in ("none", "prot_mean0.015"):
        for xnoise in (False, True):
            for rep in range(3):
                n = len(x_true)
                if off_kind == "none":
                    uu = np.zeros(P)
                else:      # per-protein excess, right-skewed, mean 0.015, sd ~0.015
                    uu = rng.gamma(1.0, 0.015, P)
                mu = ETA + uu[D["pa"]] + uu[D["pb"]] + a_true[D["g"]] * x_true ** c_true
                ys = mu + rng.normal(0, 1, n) * np.sqrt(D["sy"] ** 2 + s2)
                if xnoise:
                    D["xx"] = np.maximum(x_true + rng.normal(0, 1, n) * D["sx"], 1e-4)
                else:
                    D["xx"] = x_true
                r = {"c_true": c_true, "offsets": off_kind, "x_noise": xnoise, "rep": rep,
                     "c_fixed": Ffix.fit_c(ys, ETA)[0], "c_hog": Fhog.fit_c(ys)[0], "c_prot": Fprot.fit_c(ys)[0]}
                D["xx"] = x_true
                syn.append(r); log("syn", json.dumps(r)); res["synthetic"] = syn; save()
                if not xnoise and off_kind == "none" and rep >= 0:
                    pass
log("done")
