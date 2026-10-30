"""E4 bias diagnostics (agents/bias_fix_report.md §B): null / age-effect worlds on the euk20 TEST panel with the
real tune species sets (as scripts/tuning/v3_tuned.py part 2), tuned noise, branch engine, cells, MixedLM point
estimate (fast REML; statsmodels-equivalent). Switches turn candidate causes on/off.

Sim switches : kappa_mult (noise scale), tip_sd (per-protein tip noise; tuned value if None), nu_err (nu used by
               the engine minus the true nu), age_beta
Est switches : min_aa (merge threshold), cell_rule = aa (plan) | info (AA-only min-information rule, see
               build_cells), info_T, ycorr = none | delta (2nd-order bias-corrected log ratio), zcov = on | off,
               drop_cens
One JSON line per (config, rep) -> --out (resumable).
"""
from __future__ import annotations

import argparse
import json
import time
import warnings
from pathlib import Path

import bias_common as bc
import numpy as np

warnings.simplefilter("ignore")
import tuning_common as tc  # noqa: E402
from strucevo.stats.branch import fit_branches  # noqa: E402
from strucevo.stats.endpoints_branch import edge_view, build_cells, cell_design  # noqa: E402
from strucevo.stats.reml import fast_mixed  # noqa: E402
from strucevo.stats.synth import World, true_pair_distances, _panel_edge_lengths_for_hog  # noqa: E402
from strucevo.stats.trees import prune  # noqa: E402

DEF = dict(kappa_mult=1.0, tip_sd=None, nu_err=0.0, age_beta=0.0, nb=400, min_aa=0.05, cell_rule="aa",
           info_T=None, ycorr="none", zcov="on", drop_cens=False, struct="pair", edge_var=False)


def observe_edge(w, tuned, rng, kappa_mult=1.0, L_range=(150, 500)):
    """Edge-level noise world: realised substitutions per HOG edge (Var = kappa^2 SE(b_e; L_h)^2, AA/3Di correlated
    rho) shared by all pairs through the edge, + the non-additive pair remainder, + per-protein tip noise; Fisher SE
    columns at the observed distances (ML-engine convention), uncalibrated (calibrated below as in the driver)."""
    from strucevo.stats.synth import se_of
    k = tuned["summary"]["kappa_tdi"] * kappa_mult
    rho = tuned["summary"]["rho_tdi"]; nu = tuned["nu_tdi"]; tsd = tuned["tip_sd"]
    cols = {c: [] for c in ("hog", "species_a", "species_b", "aa_ml", "aa_ml_se", "tdi_ml", "tdi_ml_se", "aa_ml_seq")}
    for hid in w.hogs:
        ia, ib, Daa, D3, be_aa, be_3, ht = w.pairs_true[hid]
        S = w.species[hid]; n = len(ia)
        Xe = (ht.member[ia] ^ ht.member[ib]).astype(float)
        Lh = float(rng.integers(*L_range))
        va = (k * se_of("LG", be_aa, Lh)) ** 2; vt = (k * se_of("Q3Di", be_3, Lh)) ** 2
        za = rng.standard_normal(len(be_aa)); zt = rho * za + np.sqrt(1 - rho ** 2) * rng.standard_normal(len(be_aa))
        u = nu + tsd * rng.standard_normal(len(S))
        yt = D3 + u[ia] + u[ib]
        sa_t = k * se_of("LG", Daa, Lh); st_t = k * se_of("Q3Di", yt, Lh)
        z1 = rng.standard_normal(n); z2 = rho * z1 + np.sqrt(1 - rho ** 2) * rng.standard_normal(n)
        x = Daa + Xe @ (np.sqrt(va) * za) + np.sqrt(np.maximum(0, sa_t ** 2 - Xe @ va)) * z1
        y = yt + Xe @ (np.sqrt(vt) * zt) + np.sqrt(np.maximum(0, st_t ** 2 - Xe @ vt)) * z2
        cols["hog"].append(np.full(n, hid)); cols["species_a"].append(S[ia]); cols["species_b"].append(S[ib])
        cols["aa_ml"].append(x); cols["tdi_ml"].append(y)
        cols["aa_ml_se"].append(k * se_of("LG", np.maximum(x, 0), Lh)); cols["tdi_ml_se"].append(k * se_of("Q3Di", np.maximum(y, 0), Lh))
        cols["aa_ml_seq"].append(x); cols.setdefault("len_a", []).append(np.full(n, Lh))
    t = {c: np.concatenate(v) for c, v in cols.items()}
    n = len(t["hog"])
    t.update({"split": np.full(n, "tune"), "status": np.full(n, "ok"), "aa_sat": np.zeros(n, bool),
              "tdi_sat": np.zeros(n, bool), "one_to_one": np.ones(n, bool)})
    return t


def make_world(panel, Sets, nb, age_beta, rng):
    N = panel.n_nodes
    lmult = np.exp(0.3 * rng.standard_normal(N) - 0.045)
    pick = rng.integers(0, len(Sets), nb)
    hogs, species, baa, rh = [], {}, {}, {}
    for j, k in enumerate(pick):
        S = np.sort(Sets[k])
        ht0 = prune(panel, S)
        rate = rng.uniform(0.5, 1.9) / (2 * panel.age[ht0.root])
        _, b = _panel_edge_lengths_for_hog(panel, S, rate, lmult, 0.2, rng)
        hid = f"W{j:05d}"
        hogs.append(hid); species[hid] = S; baa[hid] = b; rh[hid] = 0.35 * np.exp(0.3 * rng.standard_normal())
    w = World(panel=panel, hogs=hogs, species=species, b_aa=baa, r_h=rh, lineage_mult=lmult,
              spiked=np.ones(N), age_beta=age_beta)
    for hid in hogs:
        w.pairs_true[hid] = true_pair_distances(w, hid)
    return w


def one(cfg, i, panel, tuned, Sets):
    rng = np.random.default_rng(np.random.SeedSequence(20261030, spawn_key=(2101, i)))
    w = make_world(panel, Sets, cfg["nb"], cfg["age_beta"], rng)
    tu = dict(tuned); tu["summary"] = dict(tuned["summary"])
    tu["summary"]["kappa_tdi"] = tuned["summary"]["kappa_tdi"] * cfg["kappa_mult"]
    if cfg["tip_sd"] is not None:
        tu["tip_sd"] = cfg["tip_sd"]
    t = (observe_edge(w, tu, rng, kappa_mult=cfg["kappa_mult"]) if cfg["struct"] == "edge"
         else tc.observe_world_tuned(w, tu, rng))
    nu = tuned["nu_tdi"] + cfg["nu_err"]
    bf = fit_branches(t, panel, metrics=["aa_ml", "tdi_ml"], nu={"tdi_ml": nu}, se_method="model", misfit=False,
                      display_nnls=False)
    ev = edge_view(bf, panel)
    if cfg["edge_var"]:
        from strucevo.stats.synth import se_of
        hog = t["hog"].astype(str); L = t["len_a"].astype(float)
        Lh = {h: float(np.median(L[hog == h])) for h in np.unique(hog)}
        Le = np.array([Lh[h] for h in ev.hog_ids])[ev.hog_code]
        rbar = float(np.nansum(ev.b_t) / np.nansum(ev.b_aa))
        kt = tuned["summary"]["kappa_tdi"]
        ev.se_t = np.sqrt(np.nan_to_num(ev.se_t) ** 2 + (kt * se_of("Q3Di", rbar * np.maximum(ev.b_aa, 0), Le)) ** 2)
    kw = {}
    if cfg["cell_rule"] != "aa":
        kw = {"rule": cfg["cell_rule"], "info_T": cfg["info_T"], "tip_var": tuned["tip_sd"] ** 2}
    elif cfg["ycorr"] != "none":
        kw = {"tip_var": tuned["tip_sd"] ** 2}
    cells = build_cells(ev, panel.age, min_aa=cfg["min_aa"], **kw)
    if cfg["drop_cens"]:
        m = ~cells["censored"]; cells = {k: (v[m] if isinstance(v, np.ndarray) else v) for k, v in cells.items()}
    y, x, z, _, _ = cell_design(cells, ycorr=cfg["ycorr"]) if cfg["ycorr"] != "none" else cell_design(cells)
    if cfg["zcov"] == "off":
        z = np.zeros_like(z)
        f = fast_mixed(y, x, z + 1e-9 * rng.standard_normal(len(z)), cells["hog"])
    else:
        f = fast_mixed(y, x, z, cells["hog"])
    young = cells["bin"] <= 1
    return {"beta": float(f["beta"]), "gamma": float(f.get("gamma", np.nan)), "n_cells": int(len(y)),
            "cens_frac": float(cells["censored"].mean()), "converged": bool(f.get("converged", True)),
            "y_young_minus_old": float(y[young].mean() - y[~young].mean()) if young.any() and (~young).any() else None}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", action="append", required=True)
    ap.add_argument("--reps", type=int, default=16)
    ap.add_argument("--out", default=str(bc.TEST_ROOT / "results" / "tuning" / "bias_e4_diag.jsonl"))
    a = ap.parse_args()
    D = bc.load_designs()
    Sets = [d["S"] for d in D["designs"]]
    panel, _ = tc.load_panel(bc.TEST_ROOT)
    tuned = tc.load_tuned(bc.TEST_ROOT, "tdiaa")
    done = set()
    outp = Path(a.out)
    if outp.exists():
        for line in outp.read_text().splitlines():
            r = json.loads(line); done.add((r["config"], r["rep"]))
    for spec in a.config:
        name, js = spec.split(":", 1)
        cfg = {**DEF, **json.loads(js)}
        t0 = time.time()
        for i in range(a.reps):
            if (name, i) in done:
                continue
            r = one(cfg, i, panel, tuned, Sets)
            r.update({"config": name, "rep": i, "cfg": cfg})
            with open(outp, "a") as fh:
                fh.write(json.dumps(r) + "\n")
        print(f"{name}: {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
