"""§3.2 B_E4 = the largest of {1000, 500, 200} keeping the E4 bootstrap <= 1 h.

Times run_E4 (fast REML in replicates, D6) per bootstrap replicate
  (a) on the real tune cells (branch engine on split == 'tune', --aln, tau_misfit from tau_misfit.json if present);
  (b) at production scale: a synthetic world of --prod-hogs HOGs (default 5600 = the confirm set; 8-30 species) on
      the REAL panel tree with the tuned noise (as scripts/verify/bench_stats.py).
Per-replicate cost = (t(B = b2) - t(B = b1)) / (b2 - b1) with sensitivities off. The full E4 bootstrap per alignment
also runs 3 bootstrap sensitivities with max(20, B/5) replicates each (run_E4), so
  total(B) = (B + 3 * max(20, B // 5)) * per_rep.
B_E4 is chosen on the production-scale estimate (the binding case). Threads used are recorded (this node: 2).
Writes ROOT/results/tuning/b_e4_timing.json.
"""
import argparse
import time

import numpy as np

import tuning_common as tc
from strucevo.stats.branch import fit_branches
from strucevo.stats.endpoints_branch import run_E4
from strucevo.stats.synth import make_world

ap = tc.add_common_args(argparse.ArgumentParser())
ap.add_argument("--b1", type=int, default=2)
ap.add_argument("--b2", type=int, default=12)
ap.add_argument("--prod-hogs", type=int, default=5600)
ap.add_argument("--prod-b2", type=int, default=7)
ap.add_argument("--budget-s", type=float, default=3600.0)
ap.add_argument("--skip-prod", action="store_true")
a = ap.parse_args()
t0 = time.time()
root = a.root
tuned = tc.load_tuned(root, a.aln)
tau = tc.tau_from_fragment(root)


def per_rep(bf, panel, b1, b2):
    t = time.time(); r1 = run_E4(bf, panel, tau_misfit=tau, B=b1, sensitivities=False, tip_var=tuned["tip_sd"] ** 2); t1 = time.time() - t
    t = time.time(); r2 = run_E4(bf, panel, tau_misfit=tau, B=b2, sensitivities=False, tip_var=tuned["tip_sd"] ** 2); t2 = time.time() - t
    pr = max((t2 - t1) / (b2 - b1), 1e-6)
    return {"per_rep_s": pr, "point_s": max(t1 - b1 * pr, 0.0), "t_B1_s": t1, "t_B2_s": t2, "B1": b1, "B2": b2,
            "n_hogs": r2["n_hogs"], "n_cells": r2["n_units"], "beta": r2["estimate"], "se_boot_smallB": r2["se_boot"],
            "flags": r2["flags"]}


def total(B, pr, point):
    return point + (B + 3 * max(20, B // 5)) * pr


def choose(pr, point):
    for B in (1000, 500, 200):
        if total(B, pr, point) <= a.budget_s:
            return B
    return 200


rec = tc.header(root, "b_e4_timing", a)
rec["tau_misfit_used"] = tau
bf, panel, tab, info = tc.fit_tune_branches(root, a.aln, tuned, tau_misfit=tau, metrics=["aa_ml", "tdi_ml"])
tune = per_rep(bf, panel, a.b1, a.b2)
tune["fit_info"] = info
tune["projected_total_s"] = {str(B): total(B, tune["per_rep_s"], tune["point_s"]) for B in (1000, 500, 200)}
rec["tune"] = tune
print("tune:", tune["per_rep_s"], "s/rep", tune["n_cells"], "cells")
if not a.skip_prod:
    t = time.time()
    w = make_world(panel, a.prod_hogs, rng=tc.tuning_rng("b_e4_world", 0), min_sp=8, max_sp=30)
    tabp = tc.observe_world_tuned(w, tuned, tc.tuning_rng("b_e4_world", 1))
    t_world = time.time() - t
    t = time.time()
    bfp = fit_branches(tabp, panel, metrics=["aa_ml", "tdi_ml"], nu={"tdi_ml": tuned["nu_tdi"]}, se_method="model")
    t_fit = time.time() - t
    prod = per_rep(bfp, panel, a.b1, a.prod_b2)
    prod.update({"n_pairs": int(len(tabp["hog"])), "world_s": t_world, "branch_fit_s": t_fit,
                 "panel": "real panel tree (ROOT/data/interim/timetree), ages NA-filled"})
    prod["projected_total_s"] = {str(B): total(B, prod["per_rep_s"], prod["point_s"]) for B in (1000, 500, 200)}
    rec["production_scale"] = prod
    print("prod:", prod["per_rep_s"], "s/rep", prod["n_cells"], "cells")
    basis = prod
    rec["basis"] = "production_scale"
else:
    basis = tune
    rec["basis"] = "tune (production scale skipped)"
rec["B_E4"] = choose(basis["per_rep_s"], basis["point_s"])
rec["rule"] = "largest of {1000,500,200} with point + (B + 3*max(20,B//5)) * per_rep <= budget_s (per alignment)"
rec["note"] = ("timed with OMP/torch threads = 2 on a shared node; production uses 48 threads (fast REML is "
               "numpy/BLAS-bound, so 48 threads can only be faster); the main + 2 sensitivity alignments run "
               "sequentially -> 3x per endpoint run")
rec["runtime_s"] = time.time() - t0
tc.dump(rec, tc.out_path(root, "b_e4_timing"))
print("B_E4 =", rec["B_E4"], "runtime", rec["runtime_s"])
