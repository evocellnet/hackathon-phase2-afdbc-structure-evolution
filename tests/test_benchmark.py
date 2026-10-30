"""Throughput benchmark: ML distance for 1,000,000 count matrices on GPU (LG+G4, alpha=1).

Count matrices are realistic (simulated under LG+G4 at d ~ U(0.01, 4), length ~ U(50, 800)) rather than
uniform random noise, so Newton iterations reflect real use. Result -> data/interim/tests/results/benchmark.json
"""
import json
import os
import time

import numpy as np
import pytest
import torch
from conftest import GLOBAL_SEED, RESULTS

from strucevo.distances import load_model, ml_distance, simulate_counts

pytestmark = pytest.mark.slow


@pytest.mark.skipif(not torch.cuda.is_available(), reason="needs GPU")
def test_throughput_1M():
    m = load_model("LG")
    rng = np.random.default_rng(GLOBAL_SEED)
    N = 1_000_000
    d = rng.uniform(0.01, 4.0, N)
    L = rng.integers(50, 800, N)
    t0 = time.time()
    C = simulate_counts(m, d, L, 1.0, rng=rng).astype(np.int32)
    t_sim = time.time() - t0
    Ct = torch.as_tensor(C)
    ml_distance(Ct[:10_000], m, alpha=1.0, device="cuda")      # warm-up
    out = {}
    for chunk in (100_000, 200_000, 400_000):
        torch.cuda.synchronize()
        t0 = time.time()
        r = ml_distance(Ct, m, alpha=1.0, device="cuda", chunk=chunk)
        torch.cuda.synchronize()
        dt = time.time() - t0
        out[chunk] = dict(seconds=dt, pairs_per_s=N / dt)
    rel = r.d / d - 1
    res = dict(n_pairs=N, gpu=torch.cuda.get_device_name(0), dtype="float64", sim_seconds_cpu=t_sim,
               by_chunk=out, best_pairs_per_s=max(v["pairs_per_s"] for v in out.values()),
               frac_converged=float(r.converged.mean()), frac_saturated=float(r.saturated.mean()),
               median_abs_rel_err=float(np.median(np.abs(rel))))
    # CPU reference on a subset (48 threads)
    t0 = time.time()
    ml_distance(C[:20_000], m, alpha=1.0, device="cpu")
    res["cpu_pairs_per_s"] = 20_000 / (time.time() - t0)
    json.dump(res, open(os.path.join(RESULTS, "benchmark.json"), "w"), indent=1)
    print(json.dumps(res, indent=1))
    assert res["frac_converged"] > 0.999
