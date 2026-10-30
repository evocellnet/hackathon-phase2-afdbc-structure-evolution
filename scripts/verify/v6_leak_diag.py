"""D30 diagnostic: E4 point estimates under D13 edge-tau permutations in the V6 world (true age effect), for
several cell rules (agents/bias_fix_report.md). Point estimates only (fast REML)."""
import json, sys, warnings
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
warnings.simplefilter("ignore")
from strucevo.stats.common import stage_rng
from strucevo.stats.synth import random_panel, make_world, observe_gauss
from strucevo.stats.branch import fit_branches
from strucevo.stats.endpoints_branch import edge_view, build_cells, cell_design
from strucevo.stats.reml import fast_mixed

TIP = 0.01
panel = random_panel(150, rng=stage_rng("v4_sim", 10))
w = make_world(panel, 600, rng=stage_rng("v4_sim", 11), age_beta=-0.15)
tab = observe_gauss(w, rng=stage_rng("v4_sim", 12), rho=0.3, nu_tdi=0.02, tip_sd=TIP)
bf = fit_branches(tab, panel, metrics=["aa_ml", "tdi_ml"], nu={"tdi_ml": 0.02}, se_method="model", misfit=False)
ev = edge_view(bf, panel)
nE = len(ev.child)
rules = {"plan": {}, "info4": {"rule": "info", "info_T": 4, "tip_var": TIP ** 2},
         "info4_tipall": {"rule": "info", "info_T": 4, "tip_var": TIP ** 2, "tip_all": True},
         "info3": {"rule": "info", "info_T": 3, "tip_var": TIP ** 2}}
out = {}
for name, kw in rules.items():
    c = build_cells(ev, panel.age, **kw); y, x, z, _, _ = cell_design(c); b0 = fast_mixed(y, x, z, c["hog"])["beta"]
    bs = []
    for k in range(20):
        rng = stage_rng("age_perm", k)
        c = build_cells(ev, panel.age, tau_perm=rng.permutation(nE), **kw)
        y, x, z, _, _ = cell_design(c); bs.append(fast_mixed(y, x, z, c["hog"])["beta"])
    out[name] = {"true_age_beta": float(b0), "perm_mean": float(np.mean(bs)), "perm_sd": float(np.std(bs, ddof=1)),
                 "n_cells_true": int(len(y))}
    print(name, out[name], flush=True)
