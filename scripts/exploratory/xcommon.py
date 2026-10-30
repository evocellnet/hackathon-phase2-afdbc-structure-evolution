"""Shared helpers for the exploratory analyses X1, X2, X7 (ANALYSIS_PLAN §6; unadjusted 95% CIs).

Conventions
- Only split == 'tune' is read (strucevo.stats.common.load_pairs raises for 'confirm'; guard_split re-checked).
- Primary alignment = tdiaa (DEVIATIONS D22); seqonly / tmalign are sensitivity runs (--aln).
- Seeds: np.random.SeedSequence(20261030, spawn_key=(EXPLORATORY_SPAWN[name], ...)). config/seeds.yaml is
  frozen ("never change") and has no exploratory stage, so exploratory streams use spawn keys >= 1000, which
  cannot collide with the endpoint stages (0..15).
- Tuned values (eta, kappa, nu) are read from ROOT/config/tuning.yaml; a null value falls back to the stored
  pair-table value (alpha = 1 placeholder, kappa = 1, eta = 0) and the fallback is recorded in the JSON.
- Licence: ages are used in computation only; JSON outputs hold binned summaries and fit parameters.
"""
from __future__ import annotations

import json
import os
import re
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "POLARS_MAX_THREADS"):
    os.environ.setdefault(_v, "2")        # node constraint: <= 2 CPU threads
import numpy as np
import polars as pl

R_DEFAULT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(R_DEFAULT / "src"))

from strucevo.stats.common import load_pairs, guard_split, md5_file, GLOBAL_SEED  # noqa: E402
from strucevo.stats.trees import PanelTree  # noqa: E402

EXPLORATORY_SPAWN = {"X1": 1001, "X2": 1002, "X7": 1007}
AGE_EDGES = [0.0, 30.0, 100.0, 300.0, 1000.0, np.inf]           # ANALYSIS_PLAN §0 (Mya)
AGE_LABELS = ["1 [0,30)", "2 [30,100)", "3 [100,300)", "4 [300,1000)", "5 [1000,inf)"]
# Okabe-Ito (colour-blind safe), fixed order
OI = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000", "#999999"]

# metric -> spec. y = sign*col + offset; se: 'col:<name>' | 'binom' (sqrt(p(1-p)/n_core), floor 1/n_core) | None
#   eta_key: key under tuning.yaml eta.<aln> (None => eta = 0: metric has no floor estimate in the plan)
#   excl: hogs.parquet exclusion family (excl_<fam>_<aln>)
METRICS = {
    "tm_short":    dict(col="tm_short", transform="1-", label="1 − TM-score (normalised by shorter protein)",
                        unit="", se=None, sat=None, eta_key=None, excl="lddt", additive=False),
    "rmsd_core":   dict(col="rmsd_core", transform=None, label="core RMSD", unit="Å", se=None, sat=None,
                        eta_key=None, excl="lddt", additive=False),
    "aa_p":        dict(col="aa_p", transform=None, label="AA p-distance (core columns)", unit="subst./site",
                        se="binom", sat=None, eta_key=None, excl="aa", additive=False),
    "tdi_p":       dict(col="tdi_p", transform=None, label="3Di p-distance (core columns)", unit="subst./site",
                        se="binom", sat=None, eta_key=None, excl="tdi", additive=False),
    "lddt_core":   dict(col="lddt_core", transform="1-", label="1 − LDDT (core residues)", unit="",
                        se=None, sat=None, eta_key=None, excl="lddt", additive=False),
    "lddt_global": dict(col="lddt_global", transform="1-", label="1 − LDDT global (raw 1 − S)", unit="",
                        se=None, sat=None, eta_key="lddt_raw", excl="lddt", additive=False),
    "tdi_ml":      dict(col="tdi_ml", transform=None, label="3Di ML distance (Q.3Di.AF+Γ4)", unit="subst./site",
                        se="col:tdi_ml_se", sat="tdi_sat", eta_key="tdi_ml", excl="tdi", additive=True),
    "lddt_D":      dict(col="lddt_D", transform=None, label="LDDT distance D = −ln((S−S0)/(1−S0))", unit="",
                        se="col:lddt_D_se", sat="lddt_sat", eta_key="lddt_D", excl="lddt", additive=True),
}
XVARS = {"aa_ml": dict(col="aa_ml", se="aa_ml_se", sat="aa_sat",
                       label="AA ML distance, core columns (LG+Γ4) [subst./site]"),
         "aa_ml_seq": dict(col="aa_ml_seq", se=None, sat=None,
                           label="AA ML distance, seqonly alignment (aa_ml_seq) [subst./site]")}


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def rng(name, *extra):
    return np.random.default_rng(np.random.SeedSequence(GLOBAL_SEED, spawn_key=(EXPLORATORY_SPAWN[name], *extra)))


def multinomial(r: np.random.Generator, n: int, B: int) -> np.ndarray:
    """(B, n) multinomial counts with n draws (HOG-cluster bootstrap weights)."""
    return r.multinomial(n, np.full(n, 1.0 / n), size=B).astype(float)


# ---------------------------------------------------------------------------------------------- config
def load_tuning(root: Path, aln: str) -> dict:
    """Tuned values for this alignment with recorded fallbacks."""
    out = {"eta": {}, "kappa": {}, "fallbacks": [], "tuning_md5": None}
    p = Path(root) / "config" / "tuning.yaml"
    y = {}
    if p.exists():
        import yaml
        y = yaml.safe_load(p.read_text()) or {}
        out["tuning_md5"] = md5_file(p)
    else:
        out["fallbacks"].append(f"{p} missing")
    eta = ((y.get("eta") or {}).get(aln) or {})
    for k in ("tdi_ml", "lddt_D", "lddt_raw"):
        v = eta.get(k)
        out["eta"][k] = float(v) if v is not None else 0.0
        if v is None:
            out["fallbacks"].append(f"eta.{aln}.{k}=null -> 0")
    al = y.get("alpha") or {}
    out["alpha"] = {k: al.get(k) for k in ("aa", "tdi")}
    # alpha actually used = what the pair tables were computed with (shard .done sentinels)
    used = set()
    for f in sorted((Path(root) / "data" / "interim" / "pairs" / aln).glob("*.done")):
        try:
            a_ = json.loads(f.read_text()).get("alpha") or {}
            used.add((a_.get("aa"), a_.get("tdi")))
        except Exception:
            pass
    out["alpha_in_pair_tables"] = sorted(used, key=str)
    if any(v is None for v in out["alpha"].values()):
        out["fallbacks"].append("alpha=null -> stored pair-table ML distances (alpha = 1 placeholder)")
    elif used and any((u[0] != out["alpha"]["aa"]) or (u[1] != out["alpha"]["tdi"]) for u in used):
        out["fallbacks"].append(f"tuned alpha {out['alpha']} not yet applied to pair tables (tables at {sorted(used, key=str)}); "
                                "stored table ML distances used (recompute_ml.py pending) -- eta at tuned alpha "
                                "may not match these tables")
    tm = y.get("tau_misfit")
    out["tau_misfit"] = float(tm) if tm is not None else None
    if tm is None:
        out["fallbacks"].append("tau_misfit=null -> misfit exclusion not applied")
    kp = ((y.get("kappa") or {}).get(aln) or {})
    for k in ("aa_ml", "tdi_ml"):
        edges, vals = kp.get("bin_edges"), kp.get(k)
        if edges is None or vals is None:
            out["kappa"][k] = None
            out["fallbacks"].append(f"kappa.{aln}.{k}=null -> 1")
        else:
            out["kappa"][k] = (list(map(float, edges)), list(map(float, vals)))
    return out


def kappa_fn(spec):
    """kappa(d) step function from (bin_edges, values) or 1."""
    if spec is None:
        return lambda d: np.ones_like(np.asarray(d, float))
    edges, vals = np.asarray(spec[0]), np.asarray(spec[1])
    return lambda d: vals[np.clip(np.searchsorted(edges, np.asarray(d, float), side="right") - 1, 0, len(vals) - 1)]


# ---------------------------------------------------------------------------------------------- data
def load_tune_pairs(root: Path, aln: str) -> pl.DataFrame:
    df = load_pairs(aln, "tune", root=root)
    guard_split(df["split"].to_numpy())
    assert set(df["split"].unique().to_list()) == {"tune"}
    return df


def primary_1to1(df: pl.DataFrame) -> pl.DataFrame:
    return df.filter((pl.col("pair_set") == "primary") & (pl.col("one_to_one").cast(pl.Int64) == 1))


def load_hog_flags(root: Path) -> pl.DataFrame:
    """hogs.parquet restricted to split == 'tune' (lazy filter before collect; no confirm rows materialised)."""
    p = Path(root) / "data" / "interim" / "hogs.parquet"
    lf = pl.scan_parquet(p)
    cols = [c for c in lf.collect_schema().names() if c in ("hog", "split", "selected") or c.startswith("excl_")]
    h = lf.filter(pl.col("split") == "tune").select(cols).collect()
    guard_split(h["split"].to_numpy())
    return h


def excluded_hogs(hflags: pl.DataFrame, fam: str, aln: str) -> set:
    c = f"excl_{fam}_{aln}"
    if c not in hflags.columns:
        return set()
    return set(hflags.filter(pl.col(c).fill_null(False))["hog"].to_list())


def metric_values(df: pl.DataFrame, m: str, tuning: dict | None = None):
    """(y, se_or_None, sat_bool) arrays for metric m."""
    s = METRICS[m]
    y = df[s["col"]].to_numpy().astype(float)
    if s["transform"] == "1-":
        y = 1.0 - y
    se = None
    if s["se"] == "binom":
        nc = np.maximum(df["n_core"].to_numpy().astype(float), 1.0)
        se = np.maximum(np.sqrt(np.clip(y * (1 - y), 0, None) / nc), 1.0 / nc)
    elif s["se"] and s["se"].startswith("col:"):
        se = df[s["se"][4:]].to_numpy().astype(float)
        if m == "tdi_ml" and tuning is not None:
            se = se * kappa_fn(tuning["kappa"].get("tdi_ml"))(df["aa_ml"].to_numpy())
    sat = df[s["sat"]].to_numpy().astype(bool) if s["sat"] else np.zeros(len(y), bool)
    return y, se, sat


# ---------------------------------------------------------------------------------------------- panel
def load_panel(root: Path) -> PanelTree:
    """Panel tree with leaves relabelled to taxids (leaf labels are '<UPid>_<taxid>') and node_ages medians."""
    tt = Path(root) / "data" / "interim" / "timetree"
    nwk = re.sub(r"UP\d+_(\d+)", r"\1", (tt / "panel_tree.nwk").read_text().strip())
    na = pl.read_csv(tt / "node_ages.tsv", separator="\t", infer_schema_length=100000,
                     columns=["node", "median", "lo95", "hi95"])
    ages = {r[0]: (r[1], r[2], r[3]) for r in na.iter_rows()}
    ages = {k: tuple(np.nan if v is None else float(v) for v in t) for k, t in ages.items()}
    return PanelTree.from_newick(nwk, ages=ages, is_path=False)


def mrca_age(panel: PanelTree, sa: np.ndarray, sb: np.ndarray) -> np.ndarray:
    """Age (Mya) of MRCA(species_a, species_b) for every pair (computed per unique species pair)."""
    key = np.minimum(sa, sb).astype(np.int64) * 10**10 + np.maximum(sa, sb).astype(np.int64)
    uk, inv = np.unique(key, return_inverse=True)
    out = np.full(len(uk), np.nan)
    for j, k in enumerate(uk):
        a, b = int(k // 10**10), int(k % 10**10)
        if a not in panel.leaf_of_taxid or b not in panel.leaf_of_taxid:
            continue
        i = panel.leaf_of_taxid[a]
        mb = panel.clade[panel.leaf_of_taxid[b]]
        while panel.clade[i] & mb == 0:
            i = panel.parent[i]
        out[j] = panel.age[i]
    return out[inv]


def age_bin(t: np.ndarray) -> np.ndarray:
    """§0 age-bin index 0..4 (-1 for NaN)."""
    t = np.asarray(t, float)
    b = np.searchsorted(np.asarray(AGE_EDGES[1:-1]), t, side="right")
    return np.where(np.isfinite(t), b, -1)


# ---------------------------------------------------------------------------------------------- summaries
def binned_summary(x, y, edges, min_n=20):
    """Per x-bin n, median, IQR of y (bins with n < min_n reported with n only)."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    idx = np.digitize(x, edges) - 1
    rows = []
    for k in range(len(edges) - 1):
        m = idx == k
        n = int(m.sum())
        r = dict(bin_lo=float(edges[k]), bin_hi=float(edges[k + 1]), n=n)
        if n >= min_n:
            q = np.quantile(y[m], [0.25, 0.5, 0.75])
            r.update(x_median=float(np.median(x[m])), y_q25=float(q[0]), y_median=float(q[1]), y_q75=float(q[2]))
        if n:
            rows.append(r)
    return rows


def ci95(est, reps):
    reps = np.asarray(reps, float)
    reps = reps[np.isfinite(reps)]
    if reps.size < 3:
        return dict(estimate=float(est), se_boot=None, ci95_normal=None, ci95_pct=None, n_reps=int(reps.size))
    se = float(np.std(reps, ddof=1))
    return dict(estimate=float(est), se_boot=se, ci95_normal=[float(est - 1.959964 * se), float(est + 1.959964 * se)],
                ci95_pct=[float(np.quantile(reps, 0.025)), float(np.quantile(reps, 0.975))], n_reps=int(reps.size))


def _clean(o):
    if isinstance(o, dict):
        return {str(k): _clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_clean(v) for v in o]
    if isinstance(o, (np.floating, float)):
        v = float(o)
        return v if np.isfinite(v) else None
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return _clean(o.tolist())
    return o


def write_json(obj, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(_clean(obj), indent=1))
    log("wrote", path)


def out_dirs(root: Path):
    r = Path(root) / "results" / "exploratory"
    f = Path(root) / "figures" / "exploratory"
    r.mkdir(parents=True, exist_ok=True)
    f.mkdir(parents=True, exist_ok=True)
    return r, f


def mpl_style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.alpha": 0.25, "savefig.dpi": 200})
    return plt


def common_args(ap):
    ap.add_argument("--root", default=str(R_DEFAULT), help="project (or test) root")
    ap.add_argument("--aln", default="tdiaa", choices=["tdiaa", "seqonly", "tmalign"])
    ap.add_argument("--B", type=int, default=200, help="HOG-cluster bootstrap replicates")
    ap.add_argument("--min-n", type=int, default=20, help="minimum pairs per bin for a reported median")
    return ap
