"""Tests for strucevo.models: recovery of known parameters and AIC model selection on synthetic data.

Model-selection criterion (fixed before running): data from y = b x (linear through origin) or
y = a (1 - exp(-b x)) with Gaussian noise; over 100 seeds, AIC choosing between the two generating families
(lin0 vs sat) must pick the true one in >= 90% of seeds for each generator. The full 4-model confusion matrix
(lin0, lin, sat, power) is reported in data/interim/tests/results/model_selection.json (not asserted: nested
models with one extra parameter are expected to win ~16% of the time under AIC by construction).
"""
import json
import os

import numpy as np
import pytest
from conftest import GLOBAL_SEED, RESULTS

from strucevo.models import (MODEL_NAMES, best_model, compare, divergence_x, fit_model, fit_models,
                             out_of_range_error)

N, SIGMA = 100, 0.05     # 100 points (~ number of (HOG,node) medians per stratum), noise sd 0.05


def gen(kind, rng):
    x = rng.uniform(0.02, 3.0, N)
    if kind == "lin0":
        mu = 0.3 * x
    else:
        mu = 0.9 * (1 - np.exp(-0.8 * x))
    return x, mu + rng.normal(0, SIGMA, N)


def test_parameter_recovery_noiseless():
    x = np.linspace(0.01, 4, 50)
    assert fit_model("lin0", x, 0.4 * x).params["b"] == pytest.approx(0.4)
    p = fit_model("lin", x, 0.1 + 0.4 * x).params
    assert p["a"] == pytest.approx(0.1) and p["b"] == pytest.approx(0.4)
    p = fit_model("sat", x, 0.8 * (1 - np.exp(-1.3 * x))).params
    assert p["a"] == pytest.approx(0.8, rel=1e-5) and p["b"] == pytest.approx(1.3, rel=1e-5)
    p = fit_model("power", x, 0.5 * x ** 0.6).params
    assert p["a"] == pytest.approx(0.5, rel=1e-5) and p["c"] == pytest.approx(0.6, rel=1e-5)


def test_robust_starts_extreme_scales():
    x = np.linspace(1, 1000, 80)                 # very different x scale; tiny b
    y = 50 * (1 - np.exp(-0.002 * x))
    p = fit_model("sat", x, y).params
    assert p["b"] == pytest.approx(0.002, rel=1e-4)
    x = np.linspace(0.001, 0.05, 80)
    p = fit_model("sat", x, 2 * (1 - np.exp(-60 * x))).params
    assert p["b"] == pytest.approx(60, rel=1e-4)


def test_ic_values():
    rng = np.random.default_rng(GLOBAL_SEED)
    x, y = gen("lin0", rng)
    f = fit_model("lin0", x, y)
    n = len(x)
    assert f.aic == pytest.approx(n * np.log(f.sse / n) + 4)
    assert f.bic == pytest.approx(n * np.log(f.sse / n) + 2 * np.log(n))
    rows = compare(fit_models(x, y))
    assert rows[0]["delta_aic"] == 0 and sum(r["weight_aic"] for r in rows) == pytest.approx(1)


def test_out_of_range():
    rng = np.random.default_rng(GLOBAL_SEED + 1)
    x, y = gen("sat", rng)
    o = out_of_range_error(x, y, q=0.5)
    # linear-through-origin fitted on the small-x half over-predicts the saturated large-x half
    assert o["lin0"]["bias"] > 0.05
    assert o["sat"]["rmse"] < o["lin0"]["rmse"]
    assert o["lin0"]["n_train"] + o["lin0"]["n_test"] == len(x)


def test_divergence_x():
    x = np.linspace(0.01, 3, 100)
    f_lin = fit_model("lin0", x, 0.5 * x)
    f_sat = fit_model("sat", x, 1.0 * (1 - np.exp(-0.5 * x)))
    # 0.5x vs (1-exp(-0.5x)): relative gap 10% at x where 0.5x = 1.1 (1-e^{-0.5x}) -> solve numerically
    from scipy.optimize import brentq
    xs = brentq(lambda t: 0.5 * t - 1.1 * (1 - np.exp(-0.5 * t)), 0.05, 5)
    assert divergence_x(f_lin, f_sat, 0.10) == pytest.approx(xs, abs=2e-3)


def test_aic_model_selection():
    conf = {g: {m: 0 for m in MODEL_NAMES} for g in ("lin0", "sat")}
    two_way = {"lin0": 0, "sat": 0}
    bic_two_way = {"lin0": 0, "sat": 0}
    for s in range(100):
        for g in ("lin0", "sat"):
            rng = np.random.default_rng(GLOBAL_SEED + 1000 * (g == "sat") + s)
            x, y = gen(g, rng)
            conf[g][best_model(x, y)] += 1
            two_way[g] += best_model(x, y, models=("lin0", "sat")) == g
            bic_two_way[g] += best_model(x, y, models=("lin0", "sat"), criterion="bic") == g
    json.dump(dict(n_points=N, sigma=SIGMA, two_way_aic_correct=two_way, two_way_bic_correct=bic_two_way,
                   four_way_aic_confusion=conf), open(os.path.join(RESULTS, "model_selection.json"), "w"), indent=1)
    assert two_way["lin0"] >= 90 and two_way["sat"] >= 90, two_way
