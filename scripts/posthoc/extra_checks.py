#!/usr/bin/env python3
"""POST HOC (exploratory, not pre-registered): three extra robustness checks of the tip / per-protein 3Di artefact.
 1. CATH fold-class contrast on INTERNAL edges only: k_alpha vs k_beta (ratio of sums), HOG bootstrap (B=1000,
    seed 20261030), same class assignment as E3 (hogs_endpoint_cols.cath_class = driver.cath); all-edge values for
    comparison; mean HOG pLDDT by class; alpha-vs-beta within HOG-pLDDT tertiles.
 2. E6b (MDO - LDO excess) adjusted for pLDDT: per-triplet excess (production E6b definition, E1 reference curve)
    regressed at the duplication-event level on dpLDDT = pLDDT(MDO) - pLDDT(LDO), outgroup pLDDT and mean copy pLDDT;
    intercept = adjusted excess; HOG bootstrap (B=1000, seed 20261030); stratum |dpLDDT| < 2.
 3. E1 curvature in high-confidence pairs (min pLDDT of the pair >= 90 / >= 85; all for reference), production
    filters via prepare_curve_data(extra_mask=...). (a) pre-registered fixed-eta estimator: production naive two-pass
    and SIMEX (no D29 PB correction), naive HOG x species pigeonhole bootstrap; (b) HOG-baseline eta_h and per-protein
    u_i + u_j estimators (scripts/posthoc/common.Fitter); per-protein offset size.
Read-only: production code from src/strucevo/stats; no file under results/endpoints or data/ is written.
  CONFIRM_RUN=1 python scripts/posthoc/extra_checks.py -> results/posthoc/robustness_extra.json
"""
import sys, json, time, tempfile
from pathlib import Path
import numpy as np, polars as pl
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(Path(__file__).resolve().parent))
from strucevo.stats.driver import EndpointDriver
from strucevo.stats.endpoints_curve import prepare_curve_data
from strucevo.stats.endpoints_misc import run_E6
from strucevo.stats.pairwise import CurveEngine, two_pass, simex, fit_batch
from strucevo.stats.bootstrap import multinomial_weights
from strucevo.stats.common import stage_rng
from common import Fitter

SEED = 20261030
OUT = ROOT / "results" / "posthoc" / "robustness_extra.json"
t0 = time.time()
def log(*a): print(f"[{time.time()-t0:7.0f}s]", *a, flush=True)
res = {"label": "POST HOC / exploratory; not part of the pre-registered analysis", "seed": SEED, "split": "confirm"}
def save(): OUT.write_text(json.dumps(res, indent=1, default=lambda o: float(o) if np.isscalar(o) else str(o)))
def pct(a, q=(2.5, 97.5)): return [float(v) for v in np.nanpercentile(a, q)]

d = EndpointDriver(split="confirm", run_alt_alignments=False, outdir=tempfile.mkdtemp(), log=lambda *a: None)
d.load(); d.tau = d.tun.tau_misfit()
tab = d.tabs["tdiaa"]; ETA = d.tun.eta("tdiaa")["tdi_ml"]
hc = pl.read_parquet(ROOT / "data/interim/hogs_endpoint_cols_confirm.parquet")
hogs = d.hogs
hog_plddt = dict(zip(hogs["hog"].to_list(), hogs["mean_plddt"].to_list()))
cath = dict(zip(hc["hog"].to_list(), hc["cath_class"].to_list()))
e3 = json.loads((ROOT / "results/endpoints/E3.json").read_text())
e6 = json.loads((ROOT / "results/endpoints/E6.json").read_text())
e1 = json.loads((ROOT / "results/endpoints/E1.json").read_text())

# ============================================================================ 1. fold-class contrast, internal edges
E = pl.read_parquet(ROOT / "data/interim/branches/tdiaa.parquet").filter(pl.col("split") == "confirm")
F = pl.read_parquet(ROOT / "data/interim/branches/tdiaa_fit.parquet").filter(pl.col("split") == "confirm")
F2 = F.filter(pl.col("metric").is_in(["aa_ml", "tdi_ml"]))
bad = set(F2.filter((pl.col("reason").fill_null("") != "") | (pl.col("misfit_ratio").fill_nan(0) > d.tau))["hog"].to_list())
a = E.filter(pl.col("metric") == "aa_ml").select("hog", "edge_id", "is_terminal", pl.col("b").alias("b_aa"))
t = E.filter(pl.col("metric") == "tdi_ml").select("hog", "edge_id", pl.col("b").alias("b_t"))
J = a.join(t, on=["hog", "edge_id"], how="inner").filter(~pl.col("hog").is_in(list(bad)))
J = J.with_columns(pl.col("hog").replace_strict(cath, default="none").alias("cls"),
                   pl.col("hog").replace_strict(hog_plddt, default=None).cast(pl.Float64).alias("plddt"))
uh, hinv = np.unique(J["hog"].to_numpy(), return_inverse=True)
H = len(uh)
term = J["is_terminal"].to_numpy(); ba = J["b_aa"].to_numpy(); bt = J["b_t"].to_numpy()
hcls = np.array([cath.get(h, "none") for h in uh]); hpl = np.array([hog_plddt.get(h, np.nan) for h in uh], float)
def sums(mask):
    return np.bincount(hinv[mask], ba[mask], H), np.bincount(hinv[mask], bt[mask], H)
S = {"all": sums(np.ones(len(term), bool)), "internal": sums(~term), "terminal": sums(term)}
CLASSES = ["alpha", "alphabeta", "mixed", "beta"]
p1 = {"n_hogs_included": int(H), "n_edges": int(len(term)), "class_rule": "E3 (driver.cath; hogs_endpoint_cols.cath_class)"}
p1["reproduction"] = {"k_all": float(S["all"][0].sum() / S["all"][1].sum()), "E3_k_all": e3["estimate"],
                      "k_internal": float(S["internal"][0].sum() / S["internal"][1].sum()),
                      "E3_k_internal": e3["sensitivity"]["internal_edges_only"]["estimate"],
                      "k_alpha_all": None, "E3_k_alpha": e3["cath_contrast"]["k_alpha"], "E3_k_beta": e3["cath_contrast"]["k_beta"]}
byc = {}
for c in CLASSES:
    m = hcls == c
    r = {"n_hogs": int(m.sum()), "mean_hog_plddt": float(np.nanmean(hpl[m])), "median_hog_plddt": float(np.nanmedian(hpl[m]))}
    for et in ("all", "internal", "terminal"):
        A_, T_ = S[et]; r[f"k_{et}"] = float(A_[m].sum() / T_[m].sum())
    # share of 3Di branch length on terminal edges
    r["terminal_share_3di"] = float(S["terminal"][1][m].sum() / S["all"][1][m].sum())
    r["terminal_share_aa"] = float(S["terminal"][0][m].sum() / S["all"][0][m].sum())
    byc[c] = r
p1["by_class"] = byc
p1["reproduction"]["k_alpha_all"] = byc["alpha"]["k_all"]; p1["reproduction"]["k_beta_all"] = byc["beta"]["k_all"]
rng = np.random.default_rng(SEED)
W = rng.multinomial(H, np.full(H, 1 / H), size=1000).astype(float)
ma, mb = hcls == "alpha", hcls == "beta"
boot = {}
for et in ("all", "internal", "terminal"):
    A_, T_ = S[et]
    ka = (W[:, ma] @ A_[ma]) / (W[:, ma] @ T_[ma]); kb = (W[:, mb] @ A_[mb]) / (W[:, mb] @ T_[mb])
    pa, pb = byc["alpha"][f"k_{et}"], byc["beta"][f"k_{et}"]
    dl = np.log(ka) - np.log(kb)
    boot[et] = {"k_alpha": pa, "k_alpha_ci95": pct(ka), "k_beta": pb, "k_beta_ci95": pct(kb),
                "ratio_alpha_over_beta": pa / pb, "ratio_ci95": [float(np.exp(v)) for v in pct(dl)],
                "delta_log": float(np.log(pa / pb)), "delta_log_ci95": pct(dl), "delta_log_se": float(np.std(dl, ddof=1))}
p1["alpha_vs_beta_hog_boot"] = boot
p1["plddt_alpha_minus_beta"] = byc["alpha"]["mean_hog_plddt"] - byc["beta"]["mean_hog_plddt"]
# within HOG-pLDDT tertiles (tertiles over all included HOGs)
q = np.nanquantile(hpl, [1 / 3, 2 / 3]); tert = np.digitize(hpl, q)
p1["plddt_tertile_cuts"] = q.tolist()
tt = {}
for k in range(3):
    r = {}
    for c in ("alpha", "beta"):
        m = (hcls == c) & (tert == k)
        r[c] = {"n_hogs": int(m.sum()), "mean_plddt": float(np.nanmean(hpl[m]))}
        for et in ("all", "internal"):
            A_, T_ = S[et]; r[c][f"k_{et}"] = float(A_[m].sum() / T_[m].sum())
    for et in ("all", "internal"):
        r[f"ratio_{et}"] = r["alpha"][f"k_{et}"] / r["beta"][f"k_{et}"]
    tt[f"T{k + 1}"] = r
p1["within_plddt_tertile"] = tt
res["check1_cath_internal"] = p1; save()
log("check1", json.dumps(boot["internal"]), json.dumps(boot["all"]))

# ============================================================================ E1 curve (a_h) for check 2
hcx = set(hc.filter(pl.col("excl_tdi_tdiaa"))["hog"].to_list())
data, info, idx, sp = prepare_curve_data(tab, "E1", excluded_hogs=hcx)
data.eta = float(ETA)
eng = CurveEngine(data)
f = two_pass(eng)
c_E1 = e1["estimate"]
a_h = eng.a_h(c_E1, w=f["w"]).cpu().numpy()
curve = {"c": c_E1, "eta": ETA, "a_h": dict(zip(data.hog_ids.tolist(), a_h.tolist()))}
res["E1_naive_full_reproduced"] = {"c_naive": f["c"], "E1_c_naive": e1.get("c_naive"), "n_pairs": info["n_pairs_fit"]}
log("E1 naive", f["c"], e1.get("c_naive"))

# ============================================================================ 2. E6b adjusted for pLDDT
trip, pairs = d.e6_inputs("tdiaa")
rE6 = run_E6(trip, pairs, curve, B=10)
res["check2_e6b_plddt"] = {"reproduction": {"E6b": rE6["E6b"], "E6_json_E6b": e6["E6b"], "n_b": rE6["n_units"]["b"]}}
gp = {}
for df in (pairs, tab):
    for g_, p_ in zip(df["gene_a"].to_list(), df["mean_plddt_a"].to_list()): gp[g_] = p_
    for g_, p_ in zip(df["gene_b"].to_list(), df["mean_plddt_b"].to_list()): gp[g_] = p_
st = pairs["status"].to_numpy().astype(str)
look = {}
for i, (ga, gb) in enumerate(zip(pairs["gene_a"].to_list(), pairs["gene_b"].to_list())):
    if st[i] == "ok": look[(ga, gb)] = i; look[(gb, ga)] = i
XA = pairs["aa_ml_hA"].to_numpy().astype(float); XB = pairs["aa_ml_hB"].to_numpy().astype(float)
YB = pairs["tdi_ml_hB"].to_numpy().astype(float)
ah = curve["a_h"]
def resid(g1, g2, h):
    i = look.get((g1, g2))
    if i is None or not (np.isfinite(XB[i]) and np.isfinite(YB[i])): return np.nan
    return YB[i] - (ETA + ah[h] * max(XB[i], 0.0) ** c_E1)
def xA(g1, g2):
    i = look.get((g1, g2)); return XA[i] if i is not None else np.nan
rows = {}
for r in trip.iter_rows(named=True):
    h = r["hog"]
    if r["event_level"] == "leaf" or not (r["aa_pid_seq"] <= 0.95 and r["aa_ml_seq"] <= 2): continue
    if h not in ah or not np.isfinite(ah[h]): continue
    o1 = list(r["o1"] or []); o2 = list(r["o2"] or [])
    if not o1 or not o2: continue
    c1, c2 = r["copy1"], r["copy2"]
    m1 = np.nanmean([xA(c1, o) for o in o1]); m2 = np.nanmean([xA(c2, o) for o in o1])
    if not (np.isfinite(m1) and np.isfinite(m2)): continue
    mdo, ldo = (c1, c2) if m1 > m2 else (c2, c1)
    ex, op = [], []
    for o in o2:
        v = resid(mdo, o, h) - resid(ldo, o, h)
        if np.isfinite(v): ex.append(v); op.append(gp.get(o, np.nan))
    if not ex: continue
    rows.setdefault((h, r["dcg_id"]), []).append((np.mean(ex), gp.get(mdo, np.nan), gp.get(ldo, np.nan), np.nanmean(op)))
keys = list(rows)
arr = np.array([np.mean(np.array(rows[k], float), axis=0) for k in keys])   # dcg-level means (production aggregation)
ex_, pm, pld, po = arr.T
hk = np.array([k[0] for k in keys])
ok = np.isfinite(pm) & np.isfinite(pld) & np.isfinite(po)
dp = pm - pld; pc = (pm + pld) / 2
p2 = res["check2_e6b_plddt"]
p2.update({"n_dcg": int(len(ex_)), "n_dcg_with_plddt": int(ok.sum()), "E6b_own": float(ex_.mean()),
           "E6b_own_with_plddt": float(ex_[ok].mean()),
           "dplddt_mdo_minus_ldo": {"mean": float(dp[ok].mean()), "median": float(np.median(dp[ok])),
                                    "frac_negative": float((dp[ok] < 0).mean()), "sd": float(dp[ok].std())},
           "outgroup_plddt_mean": float(po[ok].mean()), "copy_plddt_mean": float(pc[ok].mean())})
ex_, dp, po, pc, hk = ex_[ok], dp[ok], po[ok], pc[ok], hk[ok]
uh2, hi2 = np.unique(hk, return_inverse=True)
po_c = po - po.mean(); pc_c = pc - pc.mean()
models = {"M0_mean": np.ones((len(ex_), 1)),
          "M1_dplddt": np.c_[np.ones(len(ex_)), dp],
          "M2_dplddt_outgroup": np.c_[np.ones(len(ex_)), dp, po_c],
          "M3_dplddt_outgroup_copymean": np.c_[np.ones(len(ex_)), dp, po_c, pc_c]}
names = {"M0_mean": ["intercept"], "M1_dplddt": ["intercept", "dplddt"],
         "M2_dplddt_outgroup": ["intercept", "dplddt", "outgroup_plddt_c"],
         "M3_dplddt_outgroup_copymean": ["intercept", "dplddt", "outgroup_plddt_c", "copy_mean_plddt_c"]}
rng = np.random.default_rng(SEED)
WB = rng.multinomial(len(uh2), np.full(len(uh2), 1 / len(uh2)), size=1000).astype(float)
def wls(X, y, w):
    Xw = X * w[:, None]; return np.linalg.solve(X.T @ Xw, Xw.T @ y)
mres = {}
for mn, X in models.items():
    b0 = wls(X, ex_, np.ones(len(ex_)))
    bb = np.array([wls(X, ex_, WB[b][hi2]) for b in range(1000)])
    mres[mn] = {nm: {"est": float(b0[j]), "ci95": pct(bb[:, j]), "ci98": pct(bb[:, j], (1, 99)), "se": float(bb[:, j].std(ddof=1))}
                for j, nm in enumerate(names[mn])}
p2["regressions_dcg_level_hog_boot"] = mres
strata = {}
for lab, m in (("abs_dplddt_lt_2", np.abs(dp) < 2), ("abs_dplddt_lt_1", np.abs(dp) < 1), ("abs_dplddt_ge_2", np.abs(dp) >= 2),
               ("mdo_lower_plddt", dp < 0), ("mdo_higher_plddt", dp > 0)):
    ww = WB[:, hi2] * m[None]
    bs = (ww @ ex_) / ww.sum(1)
    strata[lab] = {"n_dcg": int(m.sum()), "mean_excess": float(ex_[m].mean()), "ci95": pct(bs), "ci98": pct(bs, (1, 99)),
                   "mean_dplddt": float(dp[m].mean())}
p2["strata"] = strata
res["check2_e6b_plddt"] = p2; save()
log("check2", json.dumps(mres["M2_dplddt_outgroup"]), json.dumps(strata["abs_dplddt_lt_2"]))

# ============================================================================ 3. E1 curvature in high-pLDDT pairs
pmin = tab.select(pl.min_horizontal("mean_plddt_a", "mean_plddt_b")).to_series().to_numpy()
rho_pair = d._rho_pair("tdiaa", "aa_tdi", tab)
p3 = {"estimators": {"fixed_naive": "production two_pass, y = eta + a_h x^c, eta = tuned 0.0524 (no SIMEX)",
                     "fixed_simex": "production SIMEX on the same fit (B_s = 50; plan estimator WITHOUT the D29 PB bias correction)",
                     "hog": "y = eta_h + a_h x^c (verifier's HOG-baseline estimator, scripts/posthoc/common.Fitter)",
                     "prot": "y = u_i + u_j + a_h x^c (per-protein effects)"},
      "E1_json": {"estimate": e1["estimate"], "c_naive": e1.get("c_naive"), "c_simex_plan": e1.get("c_simex_plan")},
      "subsets": {}}
B_NAIVE = 200; B_FIT = 200
FL = pl.scan_parquet(str(ROOT / "data/interim/pairs_floors/tdiaa/*.parquet")).filter(
    (pl.col("pair_set") == "floor_ident") & (pl.col("status") == "ok")).collect()
fl_y = FL["tdi_ml"].to_numpy().astype(float)
fl_pm = FL.select(pl.min_horizontal("mean_plddt_a", "mean_plddt_b")).to_series().to_numpy().astype(float)
p3["floor_ident_median_by_min_plddt"] = {str(th): float(np.median(fl_y[fl_pm >= th])) for th in (0, 85, 90)}
for lab, thr in (("plddt_min_ge90", 90.0), ("plddt_min_ge85", 85.0), ("all", None)):
    mask = np.ones(tab.height, bool) if thr is None else (np.nan_to_num(pmin, nan=-1) >= thr)
    dat, inf, ix, spp = prepare_curve_data(tab, "E1", excluded_hogs=hcx, extra_mask=mask)
    dat.eta = float(ETA)
    en = CurveEngine(dat)
    fp = two_pass(en)
    rec = {"plddt_min_threshold": thr, "n_pairs": inf["n_pairs_fit"], "n_hogs": inf["n_hogs"],
           "fixed_naive": {"c": fp["c"], "boundary": bool(fp["boundary"]), "s2": float(fp["s2"])}}
    sm = simex(en, fp["w"], fp["c"], rho=np.asarray(rho_pair, float)[ix], B_s=50, rng=np.random.default_rng(SEED + 2))
    rec["fixed_simex"] = {"c": float(sm["c_simex"]), "n_clip": int(sm["n_clip"]), "simex_minus_naive": float(sm["c_simex"] - fp["c"])}
    # naive pigeonhole (HOG x species) bootstrap, production fit_batch
    ia, ib = np.unique(np.r_[spp[0], spp[1]], return_inverse=True)[1].reshape(2, -1)
    nsp = int(max(ia.max(), ib.max()) + 1)
    brng = np.random.default_rng(SEED)
    import torch
    reps = []
    for b in range(B_NAIVE):
        Wh = multinomial_weights(brng, dat.H); M = multinomial_weights(brng, nsp)
        wb = fp["w"] * torch.as_tensor(Wh[dat.g] * M[ia] * M[ib], dtype=torch.float64, device=en.dev)
        cn, _ = fit_batch(en, w=wb, window=(fp["c"] - 0.3, fp["c"] + 0.3)); reps.append(float(cn[0]))
    reps = np.array(reps)
    rec["fixed_naive"].update({"se_boot": float(reps.std(ddof=1)), "ci95": pct(reps), "B": B_NAIVE,
                               "boot": "HOG x species pigeonhole (production design)"})
    rec["fixed_simex"]["ci95_approx"] = [rec["fixed_simex"]["c"] - 1.96 * reps.std(ddof=1), rec["fixed_simex"]["c"] + 1.96 * reps.std(ddof=1)]
    rec["fixed_simex"]["ci_note"] = "naive bootstrap SE applied to the SIMEX point (SIMEX bootstrap not run)"
    log(lab, "fixed", json.dumps(rec["fixed_naive"]), json.dumps(rec["fixed_simex"]))
    # fixed eta MATCHED to the subset: median identical-sequence floor among floor pairs with min pLDDT >= thr,
    # rescaled to the tuned alpha (stored floor table median / tuned eta)
    em = ETA * float(np.median(fl_y[fl_pm >= (thr or -1)])) / float(np.median(fl_y))
    dat2, _, _, _ = prepare_curve_data(tab, "E1", excluded_hogs=hcx, extra_mask=mask)
    dat2.eta = em
    en2 = CurveEngine(dat2); fp2 = two_pass(en2)
    sm2 = simex(en2, fp2["w"], fp2["c"], rho=np.asarray(rho_pair, float)[ix], B_s=50, rng=np.random.default_rng(SEED + 2))
    rec["fixed_matched_eta"] = {"eta": em, "n_floor_pairs": int((fl_pm >= (thr or -1)).sum()), "c_naive": fp2["c"],
                                "c_simex": float(sm2["c_simex"]), "boundary": bool(fp2["boundary"])}
    xm = dat.x[dat.w0 > 0]
    rec["x_aa_ml_q10_50_90"] = np.quantile(xm, [.1, .5, .9]).tolist()
    log(lab, "matched", json.dumps(rec["fixed_matched_eta"]))
    # Fitter-based estimators on the same pairs (raw x = aa_ml)
    m0 = dat.w0 > 0; ii = ix[m0]
    ga = tab["gene_a"].to_numpy()[ii].astype(str); gb = tab["gene_b"].to_numpy()[ii].astype(str)
    ug, g = np.unique(dat.g[m0], return_inverse=True)
    ka = np.char.add(np.char.add(g.astype(str), "|"), ga); kb = np.char.add(np.char.add(g.astype(str), "|"), gb)
    up, pinv = np.unique(np.r_[ka, kb], return_inverse=True)
    pa, pb = pinv[:len(ka)], pinv[len(ka):]
    phog = np.zeros(len(up), int); phog[pa] = g; phog[pb] = g
    D = dict(x=dat.x[m0], xx=dat.x[m0], y=dat.y[m0], sx=dat.sx[m0], sy=dat.sy[m0], g=g, H=len(ug), pa=pa, pb=pb, P=len(up), phog=phog)
    w = fp["w"].cpu().numpy()[m0]                     # production pass-2 weights
    Ffix = Fitter(D, w, "fixed"); cf = Ffix.fit_c(D["y"], ETA)
    rec["fitter_fixed_check"] = {"c": cf[0], "abs_diff_vs_production_naive": abs(cf[0] - fp["c"])}
    Fh = Fitter(D, w, "hog"); ch, qh, bh = Fh.fit_c(D["y"])
    eh = Fh.offsets(ch, D["y"]); eh1 = Fh.offsets(1.0, D["y"])
    rec["hog"] = {"c": ch, "boundary": bh, "Q_c1_minus_Qmin": Fh.Q(1.0, D["y"])[0] - qh,
                  "eta_h_median_at_chat": float(np.median(eh)), "eta_h_median_at_c1": float(np.median(eh1)),
                  "eta_h_q10_q90_at_c1": np.quantile(eh1, [.1, .9]).tolist()}
    # HOG bootstrap (seed) for fixed and hog estimators (Fitter, weights x HOG multiplicity)
    hrng = np.random.default_rng(SEED + 1)
    rf, rh = [], []
    for b in range(B_FIT):
        Wh = hrng.multinomial(D["H"], np.full(D["H"], 1 / D["H"])).astype(float)
        wb = w * Wh[g]; keep = wb > 0
        Db = {k: (v[keep] if isinstance(v, np.ndarray) and len(v) == len(w) else v) for k, v in D.items()}
        _, gg = np.unique(Db["g"], return_inverse=True); Db["g"] = gg; Db["H"] = int(gg.max() + 1)
        rf.append(Fitter(Db, wb[keep], "fixed").fit_c(Db["y"], ETA)[0]); rh.append(Fitter(Db, wb[keep], "hog").fit_c(Db["y"])[0])
    rec["fixed_naive"]["hog_boot_fitter_ci95"] = pct(rf)
    rec["hog"].update({"ci95": pct(rh), "se_boot": float(np.std(rh, ddof=1)), "B": B_FIT, "boot": "HOG multinomial (Fitter)"})
    log(lab, "hog", json.dumps(rec["hog"]))
    Fp = Fitter(D, w, "prot"); cp, qp, bp = Fp.fit_c(D["y"])
    u = Fp.offsets(cp, D["y"]); u1 = Fp.offsets(1.0, D["y"])
    npair = np.bincount(np.r_[pa, pb], minlength=D["P"]); okp = npair >= 5
    rec["prot"] = {"c": cp, "boundary": bp, "Q_c1_minus_Qmin": Fp.Q(1.0, D["y"])[0] - qp, "n_proteins_ge5_pairs": int(okp.sum())}
    rec["per_protein_offset"] = {
        "u_mean_at_chat": float(u[okp].mean()), "u_median_at_chat": float(np.median(u[okp])),
        "u_mean_at_c1": float(u1[okp].mean()), "u_median_at_c1": float(np.median(u1[okp])),
        "u_q10_q90_at_c1": np.quantile(u1[okp], [.1, .9]).tolist(),
        "excess_vs_eta_half_at_c1": float(u1[okp].mean() - ETA / 2),
        "excess_vs_eta_half_at_chat": float(u[okp].mean() - ETA / 2),
        "hog_baseline_excess_eta_h_minus_eta_at_c1_median": float(np.median(eh1) - ETA)}
    rec["pair_plddt_min_mean"] = float(np.nanmean(pmin[ii]))
    p3["subsets"][lab] = rec
    res["check3_e1_high_plddt"] = p3; save()
    log(lab, "prot", json.dumps(rec["prot"]), json.dumps(rec["per_protein_offset"]))
res["runtime_s"] = time.time() - t0; save()
log("done")
