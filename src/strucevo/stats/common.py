"""Shared utilities for the statistics layer (ANALYSIS_PLAN §0, §3, §4 header, §5).

seeds, confirm-split guard (trap #8), md5/provenance, normal/percentile CIs, Holm, JSON writer.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
from pathlib import Path

import numpy as np
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[3]
GLOBAL_SEED = 20261030
ALPHA_FW = 0.05
M_ENDPOINTS = 5
ALPHA_E = ALPHA_FW / M_ENDPOINTS          # 0.01
Z98 = norm.ppf(1 - ALPHA_E)              # two-sided 98% CI half-width multiplier (= one-sided 99%)  2.3263
ENDPOINT_KEYS = ["endpoint", "estimand", "estimate", "se_boot", "ci98", "ci_pct98", "one_sided_bound",
                 "p_value", "p_holm", "decision", "n_hogs", "n_units", "n_excluded", "sensitivity",
                 "alignment_robust", "flags", "plan_md5", "tuning_md5", "commit", "input_md5"]


# ---------------------------------------------------------------------------------------------- seeds
def _seed_map():
    p = ROOT / "config" / "seeds.yaml"
    stages = {}
    if p.exists():
        txt = p.read_text()
        for line in txt.splitlines():
            if line.strip().startswith("stages:"):
                body = line.split("{", 1)[1].rsplit("}", 1)[0]
                for kv in body.split(","):
                    k, v = kv.split(":")
                    stages[k.strip()] = int(v)
    return stages


SEED_STAGES = _seed_map()


def stage_seedseq(stage: str, *extra: int) -> np.random.SeedSequence:
    """SeedSequence(20261030) spawned per stage name (config/seeds.yaml). `extra` keys sub-streams
    (e.g. replicate index) deterministically: spawn_key = (stage_idx, *extra)."""
    if stage not in SEED_STAGES:
        raise KeyError(f"unknown seed stage {stage!r}; add it to config/seeds.yaml")
    return np.random.SeedSequence(GLOBAL_SEED, spawn_key=(SEED_STAGES[stage], *extra))


def stage_rng(stage: str, *extra: int) -> np.random.Generator:
    return np.random.default_rng(stage_seedseq(stage, *extra))


# ---------------------------------------------------------------------------------------------- guard
class ConfirmSplitGuard(RuntimeError):
    pass


def confirm_allowed() -> bool:
    return os.environ.get("CONFIRM_RUN", "") == "1"


def guard_split(split_values) -> None:
    """Raise if any of `split_values` is 'confirm' and CONFIRM_RUN != 1 (ANALYSIS_PLAN §3.5, trap #8)."""
    vals = set(np.unique(np.asarray(split_values, dtype=object)).tolist()) if split_values is not None else set()
    if "confirm" in vals and not confirm_allowed():
        raise ConfirmSplitGuard("split=='confirm' data may be read only with CONFIRM_RUN=1 (ANALYSIS_PLAN trap #8)")


def load_pairs(aln: str, split: str, columns=None, root: Path | str | None = None, pair_set=None):
    """Read data/interim/pairs/{aln}/*.parquet filtered to `split` (polars DataFrame).
    Raises ConfirmSplitGuard for split=='confirm' unless CONFIRM_RUN=1 -- BEFORE touching any file."""
    import polars as pl
    if split == "confirm" and not confirm_allowed():
        raise ConfirmSplitGuard("load_pairs(split='confirm') requires CONFIRM_RUN=1")
    root = Path(root) if root else ROOT
    files = sorted((root / "data" / "interim" / "pairs" / aln).glob("*.parquet"))
    if not files:
        raise FileNotFoundError(f"no pair shards for aln={aln}")
    lf = pl.scan_parquet([str(f) for f in files]).filter(pl.col("split") == split)
    if pair_set is not None:
        ps = [pair_set] if isinstance(pair_set, str) else list(pair_set)
        lf = lf.filter(pl.col("pair_set").is_in(ps))
    if columns is not None:
        lf = lf.select(columns)
    df = lf.collect()
    guard_split(df["split"].to_numpy() if "split" in df.columns else None)
    return df


def load_branches(aln: str, split: str, root=None):
    import polars as pl
    if split == "confirm" and not confirm_allowed():
        raise ConfirmSplitGuard("load_branches(split='confirm') requires CONFIRM_RUN=1")
    root = Path(root) if root else ROOT
    df = pl.read_parquet(root / "data" / "interim" / "branches" / f"{aln}.parquet").filter(pl.col("split") == split)
    return df


# ---------------------------------------------------------------------------------------------- provenance
def md5_file(path) -> str | None:
    p = Path(path)
    if not p.exists():
        return None
    h = hashlib.md5()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def md5_arrays(*arrays) -> str:
    h = hashlib.md5()
    for a in arrays:
        a = np.ascontiguousarray(np.asarray(a))
        if a.dtype == object:
            h.update(repr(a.tolist()).encode())
        else:
            h.update(a.tobytes())
    return h.hexdigest()


def git_commit() -> str | None:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True,
                              timeout=10).stdout.strip() or None
    except Exception:
        return None


def provenance() -> dict:
    return {"plan_md5": md5_file(ROOT / "ANALYSIS_PLAN.md"),
            "tuning_md5": md5_file(ROOT / "config" / "tuning.yaml"),
            "commit": git_commit()}


# ---------------------------------------------------------------------------------------------- inference helpers
def normal_ci(est, se, z=Z98):
    return [float(est - z * se), float(est + z * se)]


def pct_ci(reps, level=0.98):
    reps = np.asarray(reps, float)
    reps = reps[np.isfinite(reps)]
    if reps.size == 0:
        return [float("nan"), float("nan")]
    a = (1 - level) / 2
    return [float(np.quantile(reps, a)), float(np.quantile(reps, 1 - a))]


def holm(pvals: dict) -> dict:
    """Holm step-down adjusted p-values (monotone, capped at 1). NaN p-values are treated as 1."""
    keys = list(pvals)
    p = np.array([1.0 if (pvals[k] is None or not np.isfinite(pvals[k])) else float(pvals[k]) for k in keys])
    m = len(p)
    order = np.argsort(p, kind="stable")
    adj = np.empty(m)
    run = 0.0
    for rank, i in enumerate(order):
        run = max(run, min(1.0, (m - rank) * p[i]))
        adj[i] = run
    return {k: float(a) for k, a in zip(keys, adj)}


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, float):
        return None if not np.isfinite(o) else o
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return _jsonable(o.tolist())
    return o


def endpoint_record(**kw) -> dict:
    """Assemble an E*.json record with exactly the §4 key set first (extra keys allowed after)."""
    rec = {k: None for k in ENDPOINT_KEYS}
    rec.update({"n_excluded": {}, "sensitivity": {}, "flags": [], "alignment_robust": None})
    prov = provenance()
    for k in ("plan_md5", "tuning_md5", "commit"):
        rec[k] = prov[k]
    rec.update(kw)
    missing = [k for k in ENDPOINT_KEYS if k not in rec]
    assert not missing, missing
    return rec


def write_json(rec: dict, path) -> Path:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(_jsonable(rec), indent=2))
    return p
