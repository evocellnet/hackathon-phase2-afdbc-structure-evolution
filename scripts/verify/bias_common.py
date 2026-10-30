"""Shared helpers for the E1/E4 bias diagnostics (agents/bias_fix_report.md).

Designs = the euk20 TEST root's split == 'tune' HOGs (species sets, primary pairs, ok mask, n_core, NNLS AA
distance profile, r_h), exactly as scripts/tuning/v3_tuned.py builds them. Never reads split == 'confirm'.
Designs are cached (pickle) under the root's results/tuning/ so diagnostics restart fast.
"""
from __future__ import annotations

import os
import pickle
import sys
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "POLARS_MAX_THREADS"):
    os.environ.setdefault(_v, "4")
PROJECT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PROJECT / "src"))
sys.path.insert(0, str(PROJECT / "scripts" / "tuning"))

import numpy as np  # noqa: E402

TEST_ROOT = PROJECT / "data" / "interim" / "itest_euk20"


def load_designs(root=TEST_ROOT, aln="tdiaa"):
    cache = Path(root) / "results" / "tuning" / f"_bias_designs_v2_{aln}.pkl"
    if cache.exists():
        return pickle.loads(cache.read_bytes())
    import tuning_common as tc
    from strucevo.stats.branch import nnls_fit
    tuned = tc.load_tuned(root, aln)
    bf, panel, tab, info = tc.fit_tune_branches(root, aln, tuned, metrics=["aa_ml", "tdi_ml"], misfit=False)
    assert set(np.unique(tab["split"].astype(str))) == {"tune"}
    ka, kt = bf.metrics.index("aa_ml"), bf.metrics.index("tdi_ml")
    ncore = tab["n_core"].astype(float)
    designs = []
    for hf in bf.hogs:
        if hf.fit["under"][ka] or hf.fit["under"][kt]:
            continue
        bn = nnls_fit(hf.X, hf.d[ka], hf.fit["w"][ka], hf.ok[ka])
        if not np.all(np.isfinite(bn)):
            continue
        sa, st = np.nansum(hf.fit["b"][ka]), np.nansum(hf.fit["b"][kt])
        r = st / sa if sa > 0 and st > 0 else np.nan
        L = ncore[hf.pair_index]
        L = np.where(np.isfinite(L) & (L >= 20), L, np.nanmedian(L))
        designs.append({"hog": hf.hog, "S": hf.ht.species, "ia": hf.ia, "ib": hf.ib, "Daa": hf.X @ bn, "r": r,
                        "ok": hf.ok[ka] & hf.ok[kt], "L": L, "X": hf.X.astype(np.int8), "b": bn,
                        "term": np.array([bool(hf.ht.is_terminal_param[k]) if hasattr(hf.ht, "is_terminal_param") else False
                                          for k in range(hf.X.shape[1])])})
    rs = np.array([d["r"] for d in designs])
    r_med = float(np.nanmedian(rs))
    for d in designs:
        d["r"] = float(np.clip(d["r"], 0.05, 3.0)) if np.isfinite(d["r"]) else r_med
    tun = {k: tuned[k] for k in ("kappa_aa", "kappa_tdi", "rho_tdi")}
    out = {"designs": designs, "tuned_summary": tuned["summary"], "nu": tuned["nu_tdi"], "tip_sd": tuned["tip_sd"],
           "kappa_aa": (tun["kappa_aa"].edges.tolist(), tun["kappa_aa"].values.tolist()),
           "kappa_tdi": (tun["kappa_tdi"].edges.tolist(), tun["kappa_tdi"].values.tolist()),
           "rho_tdi": (tun["rho_tdi"].edges.tolist(), tun["rho_tdi"].values.tolist())}
    cache.write_bytes(pickle.dumps(out))
    return out


def bin_lookup(spec, x):
    edges, vals = np.asarray(spec[0], float), np.asarray(spec[1], float)
    k = np.clip(np.searchsorted(edges, np.nan_to_num(np.asarray(x, float), nan=0.0), side="right") - 1, 0, len(vals) - 1)
    return vals[k]
