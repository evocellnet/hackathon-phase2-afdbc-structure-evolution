"""Shared helpers for the alpha / SE-calibration / floors / apply tuning scripts (tune_alpha.py,
profile_alpha_validation.py, se_calibration.py, floors_tuned.py, apply_tuning.py).

Every script takes --root (default = project root R) and resolves all paths relative to it. Only split == 'tune'
rows are read (strucevo.stats.common.load_pairs). CPU threads are capped at 2 (shared node); heavy numerics run on
the GPU through strucevo.distances.ml_distance. Fragments put their config/tuning.yaml values at the TOP level of
ROOT/results/tuning/<step>.json with the same nesting as tuning.yaml; details go under "details".
"""
from __future__ import annotations

import glob
import json
import os
import sys
import time
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "POLARS_MAX_THREADS"):
    os.environ.setdefault(_v, "2")

PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT / "src"))

import numpy as np  # noqa: E402
import polars as pl  # noqa: E402
import torch  # noqa: E402
from scipy.optimize import minimize_scalar  # noqa: E402

torch.set_num_threads(2)

from strucevo.distances import ml_distance, profile_alpha  # noqa: E402
from strucevo.pairs.ml import D_MAX, models  # noqa: E402
from strucevo.stats.common import load_pairs, md5_file, provenance  # noqa: E402

ALNS = ("tdiaa", "seqonly", "tmalign")      # tdiaa first = primary (D22)
COARSE_ALPHAS = (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.7, 2.0, 2.5, 3.0, 4.0, 6.0, 10.0, 20.0,
                 50.0, 200.0, 1000.0)   # 1000 ~ no rate heterogeneity (boundary)
LEN_BIN_LO = (0, 100, 200, 300, 500, 800)  # floors.json bins (len_a)


def add_root_arg(ap):
    ap.add_argument("--root", default=str(PROJECT), help="project root (all paths resolved relative to it)")
    return ap


def out_json(root, step) -> Path:
    p = Path(root) / "results" / "tuning" / f"{step}.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return jsonable(o.tolist())
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def write_json(obj, path):
    tmp = str(path) + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(jsonable(obj), fh, indent=1)
    os.replace(tmp, path)
    return path


def split_md5(root):
    return md5_file(Path(root) / "data" / "interim" / "split.tsv")


def read_alpha(root, alpha_aa=None, alpha_tdi=None):
    """Tuned alpha: explicit values > ROOT/config/tuning.yaml (non-null) > ROOT/results/tuning/alpha.json."""
    if alpha_aa is not None and alpha_tdi is not None:
        return {"aa": float(alpha_aa), "tdi": float(alpha_tdi)}, "cli"
    import yaml
    y = Path(root) / "config" / "tuning.yaml"
    if y.exists():
        al = (yaml.safe_load(open(y)) or {}).get("alpha") or {}
        if al.get("aa") is not None and al.get("tdi") is not None:
            return {"aa": float(al["aa"]), "tdi": float(al["tdi"])}, str(y)
    j = out_json(root, "alpha")
    if j.exists():
        al = json.load(open(j))["alpha"]
        return {"aa": float(al["aa"]), "tdi": float(al["tdi"])}, str(j)
    raise RuntimeError("no tuned alpha (give --alpha-aa/--alpha-tdi, fill tuning.yaml or run tune_alpha.py)")


def tune_shard_files(root, aln, sub="pairs"):
    """Parquet shards of data/interim/<sub>/<aln> that contain split == 'tune' rows (reads only `split`)."""
    out = []
    for f in sorted(glob.glob(str(Path(root) / "data" / "interim" / sub / aln / "shard_[0-9][0-9][0-9][0-9].parquet"))):
        if pl.scan_parquet(f).select((pl.col("split") == "tune").any()).collect().item():
            out.append(f)
    return out


def load_tune_counts(root, aln, pair_ids, keys=("aa", "tdi"), sub="pairs"):
    """Count matrices (N,20,20) for the given tune pair_ids, in the given order, from counts_shard_*.npz."""
    want = {p: i for i, p in enumerate(pair_ids)}
    out = {k: np.zeros((len(pair_ids), 20, 20), np.int64) for k in keys}
    found = np.zeros(len(pair_ids), bool)
    for f in tune_shard_files(root, aln, sub):
        z = np.load(f.replace("shard_", "counts_shard_").replace(".parquet", ".npz"))
        pid = z["pair_id"]
        rows = np.array([want.get(p, -1) for p in pid])
        sel = rows >= 0
        for k in keys:
            out[k][rows[sel]] = z[k][sel]
        found[rows[sel]] = True
    if not found.all():
        raise RuntimeError(f"{aln}: counts missing for {(~found).sum()} pairs")
    return out


def total_ll(counts, model, alpha):
    """Pooled profile log-likelihood at one alpha (per-pair d profiled; d_max as in the pipeline)."""
    _, ll = profile_alpha(counts, model, alphas=(alpha,), d_max=D_MAX)
    return float(ll[0])


def fit_alpha(counts, model, coarse=COARSE_ALPHAS, xtol=0.002):
    """Coarse grid of the pooled profile likelihood, then bounded Brent refinement on log(alpha) between the
    neighbours of the grid maximum. Returns dict(alpha, loglik, ci95 (profile, drop 1.92), grid, grid_ll).
    If the maximum is at a grid end, alpha = that end and at_grid_boundary = True."""
    t0 = time.time()
    a, ll = profile_alpha(counts, model, alphas=coarse, d_max=D_MAX)
    i = int(np.argmax(ll))
    f = lambda x: -total_ll(counts, model, float(np.exp(x)))
    boundary = i in (0, len(a) - 1)
    if boundary:
        best, bll = float(a[i]), float(ll[i])
    else:
        r = minimize_scalar(f, bounds=(np.log(a[i - 1]), np.log(a[i + 1])), method="bounded", options={"xatol": xtol})
        best, bll = float(np.exp(r.x)), float(-r.fun)
    ci = []
    for side in (-1, 1):
        x0 = np.log(best)
        step, x = 0.02, x0
        while True:
            x = x + side * step
            if -f(x) < bll - 1.92 or abs(x - x0) > 3:
                break
            step *= 1.5
        u, w = x - side * step, x
        for _ in range(20):
            mid = 0.5 * (u + w)
            if -f(mid) >= bll - 1.92:
                u = mid
            else:
                w = mid
        ci.append(float(np.exp(0.5 * (u + w))))
    return dict(alpha=best, loglik=bll, ci95=ci, at_grid_boundary=boundary, grid=a.tolist(),
                grid_ll=ll.tolist(), n_pairs=int(len(counts)), n_sites=int(np.asarray(counts).sum()),
                seconds=round(time.time() - t0, 1))


def s0_table(root, aln, min_n=200, metric="LDDT_global"):
    """C28 S0(L) table for one alignment from ROOT/results/floors/floors.json: per length bin, S0_ok (random pairs
    with status ok) if it has >= min_n pairs, else the unconditional S0 of that bin. tuning.yaml layout:
    {bin_lo, median, q95, source} (+ n)."""
    fl = json.load(open(Path(root) / "results" / "floors" / "floors.json"))
    f = fl[aln]
    lo, med, q95, srcs, ns = [], [], [], [], []
    for (b0, _), name in zip(fl["bin_edges_len_a"], fl["bins"]):
        ok = f["S0_ok"][name][metric]
        un = f["S0"][name][metric]
        use, src = (ok, "S0_ok") if ok["n"] >= min_n else (un, "S0")
        lo.append(int(b0)); med.append(float(use["median"])); q95.append(float(use["q95"]))
        srcs.append(src); ns.append(int(use["n"]))
    return {"bin_lo": lo, "median": med, "q95": q95, "source": srcs, "n": ns}


def read_s0(root, aln):
    """S0 table: ROOT/config/tuning.yaml s0[aln] if filled, else ROOT/results/tuning/floors_tuned.json, else built
    from floors.json (C28)."""
    import yaml
    y = Path(root) / "config" / "tuning.yaml"
    if y.exists():
        s = ((yaml.safe_load(open(y)) or {}).get("s0") or {}).get(aln) or {}
        if s.get("bin_lo") is not None and s.get("median") is not None:
            return s, str(y)
    j = out_json(root, "floors_tuned")
    if j.exists():
        return json.load(open(j))["s0"][aln], str(j)
    return s0_table(root, aln), "floors.json (C28 built)"


def bin_index(x, edges_lo):
    """Index of the bin with the largest lower edge <= x (clipped to the first / last bin)."""
    lo = np.asarray(edges_lo, float)
    return np.clip(np.searchsorted(lo, np.asarray(x, float), side="right") - 1, 0, len(lo) - 1)


def fill_nearest(vals):
    """Fill None/NaN entries of a 1-D list by the nearest non-missing entry (ties -> lower bin). Returns
    (filled list, filled-flag list)."""
    v = np.array([np.nan if x is None else x for x in vals], float)
    good = np.where(np.isfinite(v))[0]
    flag = ~np.isfinite(v)
    if len(good) == 0:
        return [None] * len(v), flag.tolist()
    for i in np.where(flag)[0]:
        v[i] = v[good[np.argmin(np.abs(good - i) + 1e-3 * (good > i))]]
    return v.tolist(), flag.tolist()
