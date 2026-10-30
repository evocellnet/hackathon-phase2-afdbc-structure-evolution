"""ML distances from stored count matrices, and the floor-corrected LDDT distance (ANALYSIS_PLAN §0, §2.1).

aa_ml : LG+Γ4(α_aa);  tdi_ml : Q.3Di.AF+Γ4(α_3di), loaded with AA_ORDER (trap #4). d_max = 20 (trap #3).
SE floor: se = max(se, 1/n_sites); NaN SE (boundary d = 0 or d_max) -> 1/n_sites.
α defaults to a placeholder 1.0 until §Tuning fixes it; everything here is recomputable from counts_*.npz.
"""
from __future__ import annotations

import hashlib
import os

import numpy as np

from ..distances import AA_ORDER, load_model, load_paml, ml_distance

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
Q3DI_PATH = os.path.join(ROOT, "data/raw/models/Q.3Di.AF")
Q3DI_MD5 = "8f897fffd14f247afbee44dad2fce2ce"
D_MAX = 20.0

_MODELS = {}


def models():
    if not _MODELS:
        md5 = hashlib.md5(open(Q3DI_PATH, "rb").read()).hexdigest()
        if md5 != Q3DI_MD5:
            raise RuntimeError(f"Q.3Di.AF md5 {md5} != {Q3DI_MD5}")
        _MODELS["aa"] = load_model("LG")
        _MODELS["tdi"] = load_paml(Q3DI_PATH, alphabet=AA_ORDER, name="Q.3Di.AF")
    return _MODELS


def ml_from_counts(counts: np.ndarray, alphabet: str, alpha: float = 1.0, device: str | None = None) -> dict:
    """counts (N,20,20) -> dict of arrays: p, ml, se (floored), sat, n_sites."""
    m = models()["aa" if alphabet == "aa" else "tdi"]
    counts = np.asarray(counts)
    n = counts.sum((1, 2)).astype(float)
    if len(counts) == 0:
        e = np.zeros(0)
        return {"p": e, "ml": e, "se": e, "sat": np.zeros(0, bool), "n_sites": e}
    r = ml_distance(counts, m, alpha=alpha, ncat=4, d_max=D_MAX, device=device)
    with np.errstate(invalid="ignore", divide="ignore"):
        p = 1.0 - np.trace(counts, axis1=1, axis2=2) / n
        floor = 1.0 / n
    se = np.where(np.isnan(r.se), floor, np.fmax(r.se, floor))
    d = r.d.copy()
    nos = n == 0
    d[nos] = np.nan; se[nos] = np.nan; p[nos] = np.nan
    return {"p": p, "ml": d, "se": se, "sat": np.asarray(r.saturated, bool) & ~nos, "n_sites": n}


def lddt_D(S, len_a, s0_table, se_model=None) -> dict:
    """Floor-corrected LDDT distance D = -ln((S - S0(L)) / (1 - S0(L))) (ANALYSIS_PLAN §0, §1.4).

    s0_table: mapping with keys 'bin_lo' (ascending lower edges of length bins), 'median', 'q95' (arrays) —
    the S0 table (C28: config/tuning.yaml s0.{aln}, i.e. S0_ok with fallback). Pairs with S <= q95(bin) get
    lddt_sat = True (kept, censored; never dropped); D is NaN only where it is undefined (S <= S0 median).
    se_model (optional, §2.3; config/tuning.yaml lddt_se_model.{aln}): {'len_bin_lo', 's_bin_lo', 'sd_S'
    (len x S table of the block-bootstrap SD of S)} -> lddt_D_se = sd_S(L, S) / (S - S0(L)) (delta method,
    |dD/dS| = 1/(S - S0)); NaN where D is NaN. Without se_model lddt_D_se is NaN.
    """
    S = np.asarray(S, float)
    L = np.asarray(len_a, float)
    lo = np.asarray(s0_table["bin_lo"], float)
    b = np.clip(np.searchsorted(lo, L, side="right") - 1, 0, len(lo) - 1)
    s0 = np.asarray(s0_table["median"], float)[b]
    q95 = np.asarray(s0_table["q95"], float)[b]
    sat = S <= q95
    with np.errstate(invalid="ignore", divide="ignore"):
        D = -np.log((S - s0) / (1.0 - s0))
    D = np.where((S <= s0) | np.isnan(S), np.nan, D)
    se = np.full_like(D, np.nan)
    if se_model is not None:
        llo = np.asarray(se_model["len_bin_lo"], float)
        slo = np.asarray(se_model["s_bin_lo"], float)
        tab = np.asarray(se_model["sd_S"], float)
        bl = np.clip(np.searchsorted(llo, L, side="right") - 1, 0, len(llo) - 1)
        bs = np.clip(np.searchsorted(slo, np.nan_to_num(S), side="right") - 1, 0, len(slo) - 1)
        with np.errstate(invalid="ignore", divide="ignore"):
            se = np.where(np.isfinite(D), tab[bl, bs] / (S - s0), np.nan)
    return {"lddt_D": D, "lddt_D_se": se, "lddt_sat": sat}
