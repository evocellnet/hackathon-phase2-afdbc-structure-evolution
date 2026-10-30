"""Tests for the tuning helpers: profile_alpha recovery (§3.2), lddt_D_se delta method (§2.3), C28 S0 reader."""
import importlib.util
import json
import os

import numpy as np
import pytest

from strucevo.distances import load_model, profile_alpha, simulate_counts
from strucevo.pairs.ml import lddt_D

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


@pytest.mark.parametrize("alpha", [0.5, 2.0])
def test_profile_alpha_recovers_pooled_alpha(alpha):
    m = load_model("LG")
    rng = np.random.default_rng(7)
    d = rng.uniform(0.1, 2.0, 3000)
    C = simulate_counts(m, d, rng.integers(200, 401, 3000), alpha=alpha, rng=rng)
    grid = alpha * np.array([0.6, 0.8, 0.9, 1.0, 1.1, 1.25, 1.6])
    a, ll = profile_alpha(C, m, alphas=grid, d_max=20)
    assert a[np.argmax(ll)] == pytest.approx(alpha, rel=0.11)


def test_lddt_D_se_delta_method():
    s0 = {"bin_lo": [0, 100], "median": [0.2, 0.15], "q95": [0.3, 0.25]}
    se_model = {"len_bin_lo": [0, 100], "s_bin_lo": [0.0, 0.5], "sd_S": [[0.02, 0.01], [0.015, 0.008]]}
    S = np.array([0.4, 0.8, 0.1]); L = np.array([50, 300, 300])
    r = lddt_D(S, L, s0, se_model)
    eps = 1e-6
    dD = (lddt_D(S + eps, L, s0)["lddt_D"] - lddt_D(S - eps, L, s0)["lddt_D"]) / (2 * eps)
    sd = np.array([0.02, 0.008, np.nan])
    assert np.allclose(r["lddt_D_se"][:2], (np.abs(dD) * sd)[:2])
    assert np.isnan(r["lddt_D_se"][2]) and np.isnan(r["lddt_D"][2])
    assert np.all(np.isnan(lddt_D(S, L, s0)["lddt_D_se"]))   # no model -> NaN (backward compatible)


def test_recompute_ml_reads_c28_s0_from_floors_json():
    p = os.path.join(ROOT, "results", "floors", "floors.json")
    if not os.path.exists(p):
        pytest.skip("floors.json not built")
    spec = importlib.util.spec_from_file_location("rml", os.path.join(ROOT, "scripts", "pairs", "recompute_ml.py"))
    rml = importlib.util.module_from_spec(spec); spec.loader.exec_module(rml)
    f = json.load(open(p))
    t = rml.read_s0(p, "tdiaa")
    for i, name in enumerate(f["bins"]):
        ok = f["tdiaa"]["S0_ok"][name]["LDDT_global"]
        exp = ok if ok["n"] >= 200 else f["tdiaa"]["S0"][name]["LDDT_global"]
        assert t["median"][i] == exp["median"] and t["q95"][i] == exp["q95"]
