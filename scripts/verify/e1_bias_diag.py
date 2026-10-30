"""E1 bias diagnostics (agents/bias_fix_report.md §A): switch candidate causes on/off on euk20 TEST tune designs.

Linear (or power) truth y = eta + r_h D^c + tip noise, x = D + e_x, with realistic SE profiles (se_of, n_core),
kappa(bin), rho(bin). Point estimates only (two-pass naive + SIMEX); no bootstrap.

Simulation switches  : noise=bin|scalar ; se_at=true|obs (SE columns computed at the true or the observed
                       distance, production = obs) ; Lmult (n_core multiplier: noise magnitude) ; trunc (x<=0 /
                       negative y handled as in production) ; seqproxy=indep|shared
Estimator switches   : est_kappa=scalar|bin_obs|bin_true ; est_rho=scalar|pair ; weights=sey|xmodel|flat ;
                       x_input=raw|tree (tree-smoothed AA distance) ; Bs
Output: one JSON line per (config, replicate) appended to --out (resumable: finished keys are skipped).
"""
from __future__ import annotations

import argparse
import json
import time
import warnings
from pathlib import Path

import bias_common as bc
import numpy as np
import torch

from strucevo.stats.endpoints_curve import prepare_curve_data
from strucevo.stats.pairwise import CurveEngine, two_pass, simex
from strucevo.stats.synth import se_of

warnings.simplefilter("ignore")
torch.set_num_threads(4)

DEF_SIM = dict(noise="bin", se_at="obs", Lmult=1.0, c=1.0, seqproxy="indep", hogs_scale=4.0, eta_mode="fixed",
               struct="pair")
DEF_EST = dict(est_kappa="bin_obs", est_rho="scalar", weights="sey", Bs=20, lambdas=(0.5, 1.0, 1.5, 2.0),
               extrap="quad", x_input="raw", eta_shift=0.0, clip_mode="clip", pbK=0, pb_design="x")
_PANEL = {}
RES_RATIO = [float("nan")]
PB = [None]


def tree_smooth(D, t, idx, kx_all):
    """x_hat = X b_AA (pass-2 WLS on the HOG species tree, calibrated SEs) and its model SD / hat diagonal."""
    from strucevo.stats.branch import fit_branches
    import tuning_common as tc
    if "p" not in _PANEL:
        _PANEL["p"] = tc.load_panel(bc.TEST_ROOT)[0]
    tb = dict(t)
    tb["aa_ml_se"] = t["aa_ml_se"] * kx_all
    hogsp = {}
    for h, a_, b_ in zip(t["hog"], t["species_a"], t["species_b"]):
        hogsp.setdefault(h, set()).update((int(a_), int(b_)))
    hogsp = {h: np.array(sorted(v)) for h, v in hogsp.items()}
    bf = fit_branches(tb, _PANEL["p"], metrics=["aa_ml"], hog_species=hogsp, misfit=False, display_nnls=False,
                      se_method="model")
    n = len(t["hog"])
    xh = np.full(n, np.nan); sd = np.full(n, np.nan); hii = np.full(n, np.nan)
    for hf in bf.hogs:
        if hf.fit["under"][0]:
            continue
        X = hf.X; A = hf.fit["A"][0]; V = hf.fit["V"][0]
        xh[hf.pair_index] = X @ hf.fit["b"][0]
        sd[hf.pair_index] = np.sqrt(np.einsum("ij,jk,ik->i", X, V, X))
        hii[hf.pair_index] = np.einsum("ij,ji->i", X, A)
    return xh, sd, hii


def simulate(D, rng, sim):
    ds = D["designs"]
    nH = len(ds)
    n_per = int(round(sim["hogs_scale"] * nH))
    pick = rng.integers(0, nH, n_per)
    nu, tip_sd = D["nu"], D["tip_sd"]
    S = D["tuned_summary"]
    cols = {k: [] for k in ("hog", "species_a", "species_b", "aa_ml", "tdi_ml", "aa_ml_seq", "len_a",
                            "Daa", "Dy")}
    for j, k in enumerate(pick):
        d = ds[k]
        n = len(d["ia"])
        u = nu + tip_sd * rng.standard_normal(len(d["S"]))
        y_true = d["r"] * np.maximum(d["Daa"], 0.0) ** sim["c"] + u[d["ia"]] + u[d["ib"]]
        L = d["L"] * sim["Lmult"]
        sx = se_of("LG", d["Daa"], L); sy = se_of("Q3Di", y_true, L)
        if sim["noise"] == "bin":
            kx = bc.bin_lookup(D["kappa_aa"], d["Daa"]); ky = bc.bin_lookup(D["kappa_tdi"], d["Daa"])
            rr = bc.bin_lookup(D["rho_tdi"], d["Daa"])
        else:
            kx, ky, rr = S["kappa_aa"], S["kappa_tdi"], S["rho_tdi"]
        z1 = rng.standard_normal(n); z2 = rr * z1 + np.sqrt(1 - rr ** 2) * rng.standard_normal(n)
        if sim["struct"] == "edge":
            # realised substitutions per EDGE (shared by all pairs through it) + non-additive pair remainder:
            # Var(edge e) = kappa^2 SE(b_e; Lh)^2 ; pair remainder = max(0, total pair var - sum over path)
            X = d["X"].astype(float); b = np.maximum(d["b"], 0.0); Lh = float(np.median(L))
            kxe = bc.bin_lookup(D["kappa_aa"], b); kye = bc.bin_lookup(D["kappa_tdi"], b)
            rre = bc.bin_lookup(D["rho_tdi"], b)
            ve_a = (kxe * se_of("LG", b, Lh)) ** 2 * (b > 0)
            bt = d["r"] * b
            ve_t = (kye * se_of("Q3Di", bt, Lh)) ** 2 * (b > 0)
            za = rng.standard_normal(len(b)); zt = rre * za + np.sqrt(1 - rre ** 2) * rng.standard_normal(len(b))
            ea_edge = X @ (np.sqrt(ve_a) * za); et_edge = X @ (np.sqrt(ve_t) * zt)
            va_res = np.maximum(0.0, (kx * sx) ** 2 - X @ ve_a); vt_res = np.maximum(0.0, (ky * sy) ** 2 - X @ ve_t)
            ex_ = ea_edge + np.sqrt(va_res) * z1; ey_ = et_edge + np.sqrt(vt_res) * z2
            cols.setdefault("_vres_frac", []).append(va_res / np.maximum((kx * sx) ** 2, 1e-12))
        else:
            ex_ = kx * sx * z1; ey_ = ky * sy * z2
        cols["hog"].append(np.full(n, f"{d['hog']}#{j}"))
        cols["species_a"].append(d["S"][d["ia"]]); cols["species_b"].append(d["S"][d["ib"]])
        cols["aa_ml"].append(d["Daa"] + ex_)
        cols["tdi_ml"].append(y_true + ey_)
        zs = rng.standard_normal(n) if sim["seqproxy"] == "indep" else z1
        cols["aa_ml_seq"].append(d["Daa"] + kx * sx * zs)
        cols["len_a"].append(L); cols["Daa"].append(d["Daa"]); cols["Dy"].append(y_true)
    t = {k: np.concatenate(v) for k, v in cols.items() if v}
    n = len(t["hog"])
    t.update({"split": np.full(n, "tune"), "status": np.full(n, "ok"), "aa_sat": np.zeros(n, bool),
              "tdi_sat": np.zeros(n, bool), "one_to_one": np.ones(n, bool)})
    # Fisher SE columns: at the true distance (old simulation) or at the observed one (production ML engine)
    if sim["se_at"] == "true":
        t["aa_ml_se"] = se_of("LG", t["Daa"], t["len_a"]); t["tdi_ml_se"] = se_of("Q3Di", t["Dy"], t["len_a"])
    else:
        t["aa_ml_se"] = se_of("LG", np.maximum(t["aa_ml"], 0), t["len_a"])
        t["tdi_ml_se"] = se_of("Q3Di", np.maximum(t["tdi_ml"], 0), t["len_a"])
    return t


def estimate(D, t, est, eta, rng):
    S = D["tuned_summary"]
    data, inf, idx, sp = prepare_curve_data(t, "E1")
    x = t["aa_ml"][idx]
    # calibrated SEs (kappa) -> SIMEX gets kx = ky = 1 with these columns (driver convention)
    if est["est_kappa"] == "scalar":
        kx, ky = S["kappa_aa"], S["kappa_tdi"]
    elif est["est_kappa"] == "bin_obs":
        kx, ky = bc.bin_lookup(D["kappa_aa"], x), bc.bin_lookup(D["kappa_tdi"], x)
    else:
        kx, ky = bc.bin_lookup(D["kappa_aa"], t["Daa"][idx]), bc.bin_lookup(D["kappa_tdi"], t["Daa"][idx])
    data.sx = data.sx * kx
    data.sy = data.sy * ky
    rho_scale = 1.0
    if est["x_input"] == "tree":
        kx_all = bc.bin_lookup(D["kappa_aa"], t["aa_ml"])
        xh, sdh, hii = tree_smooth(D, t, idx, kx_all)
        xh, sdh, hii = xh[idx], sdh[idx], hii[idx]
        bad = ~(np.isfinite(xh) & (xh > 0))
        data.w0 = np.where(bad, 0.0, data.w0)
        rho_scale = np.where(bad, 0.0, hii * data.sx / np.where(bad, 1.0, sdh))
        data.x = np.where(bad, 1.0, xh); data.sx = np.where(bad, 1.0, sdh)
        x = data.x
    if est["weights"] == "xmodel":
        # SE_y predicted from x only: SE_y(r_hat_h x, L) with a per-HOG naive linear slope (no y-dependence per pair)
        g = data.g
        num = np.bincount(g, (data.y - eta) * x, data.H); den = np.bincount(g, x * x, data.H)
        ah = np.maximum(num / np.maximum(den, 1e-12), 1e-3)
        sy_mod = se_of("Q3Di", eta + ah[g] * x, t["len_a"][idx]) * (ky if np.ndim(ky) else ky)
        data_sy_w = sy_mod
    elif est["weights"] == "flat":
        data_sy_w = np.full(len(x), np.median(data.sy))
    else:
        data_sy_w = None
    if est["est_rho"] == "scalar":
        rho = S["rho_tdi"]
    elif est["est_rho"] == "pair":
        rho = bc.bin_lookup(D["rho_tdi"], x)
    else:
        rho = 0.0
    rho = rho * rho_scale
    xtree = None
    if est["x_input"] == "tree2":
        from strucevo.stats.branch import fit_branches
        from strucevo.stats.endpoints_curve import tree_x
        import tuning_common as tc
        if "p" not in _PANEL:
            _PANEL["p"] = tc.load_panel(bc.TEST_ROOT)[0]
        tb = dict(t); tb["aa_ml_se"] = t["aa_ml_se"] * bc.bin_lookup(D["kappa_aa"], t["aa_ml"])
        hogsp = {}
        for h, a_, b_ in zip(t["hog"], t["species_a"], t["species_b"]):
            hogsp.setdefault(h, set()).update((int(a_), int(b_)))
        bf = fit_branches(tb, _PANEL["p"], metrics=["aa_ml"], hog_species={h: np.array(sorted(v)) for h, v in hogsp.items()},
                          misfit=False, display_nnls=False, se_method="model")
        xtree, _ = tree_x(bf, data, idx)
        num = den = 0.0
        for hf in bf.hogs:
            if hf.fit["under"][0]:
                continue
            w_ = hf.fit["w"][0]; ok_ = w_ > 0; h_ = np.einsum("ij,ji->i", hf.X, hf.fit["A"][0])
            r_ = hf.d[0] - hf.X @ hf.fit["b"][0]
            num += float((r_[ok_] ** 2).sum()); den += float((hf.se[0][ok_] ** 2 * (1 - h_[ok_])).sum())
        RES_RATIO[0] = num / den
    data.eta = eta
    eng = CurveEngine(data)
    if data_sy_w is not None:
        sy_keep = eng.sy
        eng.sy = torch.as_tensor(data_sy_w, dtype=torch.float64, device=eng.dev)
        f = two_pass(eng)
        eng.sy = sy_keep
    else:
        f = two_pass(eng)
    sm = simex(eng, f["w"], f["c"], kx=1.0, ky=1.0, rho=rho, B_s=est["Bs"], rng=rng, lambdas=tuple(est["lambdas"]),
               xtree=xtree, clip_mode=est["clip_mode"])
    if est["pbK"]:
        from strucevo.stats.pairwise import simex_pb_correction
        xd = None
        if est["pb_design"] == "tree":
            import copy
            from strucevo.stats.branch import fit_branches
            from strucevo.stats.endpoints_curve import tree_x
            import tuning_common as tc
            if "p" not in _PANEL:
                _PANEL["p"] = tc.load_panel(bc.TEST_ROOT)[0]
            tb = dict(t); tb["aa_ml_se"] = t["aa_ml_se"] * bc.bin_lookup(D["kappa_aa"], t["aa_ml"])
            hogsp = {}
            for h, a_, b_ in zip(t["hog"], t["species_a"], t["species_b"]):
                hogsp.setdefault(h, set()).update((int(a_), int(b_)))
            bfd = fit_branches(tb, _PANEL["p"], metrics=["aa_ml"], hog_species={h: np.array(sorted(v)) for h, v in hogsp.items()},
                               misfit=False, display_nnls=False, se_method="model")
            dc = copy.deepcopy(data); tree_x(bfd, dc, idx)
            xd = np.where(dc.w0 > 0, dc.x, data.x)
        pb = simex_pb_correction(eng, f, sm["c_simex"], kx=1.0, ky=1.0, rho=rho, K=est["pbK"], B_s=est["Bs"], rng=rng,
                                 x_design=xd)
        PB[0] = pb
    lam = np.array(list(sm["c_lambda"].keys())); cl = np.array(list(sm["c_lambda"].values()))
    out = {"res_ratio": RES_RATIO[0], "c_pb": PB[0]["c_corrected"] if (est["pbK"] and PB[0]) else float("nan"),
           "pb_mc_se": PB[0]["mc_se"] if (est["pbK"] and PB[0]) else float("nan"),
           "c_naive": f["c"], "c_simex": sm["c_simex"], "c_lambda": cl.tolist(), "n_clip": sm["n_clip"],
           "n_pairs": inf["n_pairs_fit"], "s2": f["s2"]}
    # alternative extrapolants from the same lambda path
    try:
        from scipy.optimize import curve_fit
        p, _ = curve_fit(lambda l, a, b, c: a + b / (c + l), lam, cl, p0=(cl[0], 0.01, 2.0), maxfev=5000)
        out["c_rational"] = float(p[0] + p[1] / (p[2] - 1.0)) if p[2] > 1.05 else float("nan")
    except Exception:
        out["c_rational"] = float("nan")
    out["c_cubic"] = float(np.polyval(np.polyfit(lam, cl, 3), -1.0)) if len(lam) >= 5 else float("nan")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", action="append", required=True,
                    help="name:json overrides of sim/est keys, e.g. 'base:{}' or 'scal:{\"noise\":\"scalar\"}'")
    ap.add_argument("--reps", type=int, default=16)
    ap.add_argument("--out", default=str(bc.TEST_ROOT / "results" / "tuning" / "bias_e1_diag.jsonl"))
    a = ap.parse_args()
    D = bc.load_designs()
    eta = 2 * D["nu"]
    done = set()
    outp = Path(a.out)
    if outp.exists():
        for line in outp.read_text().splitlines():
            r = json.loads(line); done.add((r["config"], r["rep"]))
    for spec in a.config:
        name, js = spec.split(":", 1)
        ov = json.loads(js)
        sim = {**DEF_SIM, **{k: v for k, v in ov.items() if k in DEF_SIM}}
        est = {**DEF_EST, **{k: v for k, v in ov.items() if k in DEF_EST}}
        t0 = time.time()
        for i in range(a.reps):
            if (name, i) in done:
                continue
            rng = np.random.default_rng(np.random.SeedSequence(20261030, spawn_key=(2001, i)))
            t = simulate(D, rng, sim)
            r = estimate(D, t, est, eta + est["eta_shift"], np.random.default_rng(np.random.SeedSequence(20261030, spawn_key=(2002, i))))
            r.update({"config": name, "rep": i, "sim": sim, "est": est})
            with open(outp, "a") as fh:
                fh.write(json.dumps(r) + "\n")
        print(f"{name}: {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
