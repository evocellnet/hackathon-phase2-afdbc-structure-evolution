#!/usr/bin/env python3
"""Review: E4 point beta over a fine grid of tip corrections nu (incl. C35's nu +/- 0.005 and the
mean-based floor nu = mean(floor_ident tdi_ml)/2 = 0.0313). Point estimates only.
  CONFIRM_RUN=1 python scripts/review/e4_nu_grid.py > logs/review/e4_nu_grid.json"""
import sys, json
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "scripts" / "review"))
from e4_tip_diagnostics import load_bf, with_nu, e4_point, ksum
from strucevo.stats.driver import Tuning
tun = Tuning(ROOT); tau = tun.tau_misfit(); tv = tun.tip_sd("tdiaa") ** 2; nu0 = tun.eta("tdiaa")["tdi_ml"] / 2
bf, panel = load_bf()
out = []
for nu in (0.0212, 0.0262, 0.0287, 0.0313, 0.0337, 0.0362, 0.040):
    b, n, cf, ev, _ = e4_point(with_nu(bf, nu, nu0), panel, tau, tv)
    bi, ni, cfi, _, _ = e4_point(with_nu(bf, nu, nu0), panel, tau, tv, internal_only=True)
    out.append({"nu": nu, "beta": b, "n_cells": n, "cens": cf, "k_all": ksum(ev), "beta_internal_check": bi})
print(json.dumps(out, indent=1, default=float))
