"""V11: assemble results/verify/V11.json + agents/v11_report.md from the per-endpoint V11_*.json files."""
import numpy as np
from scipy.stats import norm
from v11_common import *

L = lambda n: json.loads((OUT / n).read_text())
c1, c2, pb1, pb2 = L("V11_E1_curve.json"), L("V11_E2_curve.json"), L("V11_E1_pb.json"), L("V11_E2_pb.json")
br, e5, e6, inp = L("V11_branch.json"), L("V11_E5.json"), L("V11_E6.json"), L("V11_inputs.json")
P = {e: endpoint(e) for e in ("E1", "E2", "E3", "E4", "E5", "E6")}
TOL = 1e-6
rows = []


def row(ep, q, pv, ov, kind, tol=TOL, note=""):
    d = abs(pv - ov)
    rows.append({"endpoint": ep, "quantity": q, "project": pv, "own": ov, "absdiff": d, "kind": kind, "tol": tol,
                 "pass": bool(d <= tol), "note": note})


# E1/E2
for ep, c, pb in (("E1", c1, pb1), ("E2", c2, pb2)):
    p = P[ep]
    row(ep, "n pairs in fit", p["n_units"], c["counts_own"]["n_fit"], "exact", 0)
    row(ep, "n HOGs in fit", p["n_hogs"], c["counts_own"]["n_hogs"], "exact", 0)
    row(ep, "s2 (pass-1 excess var)", p["s2"], c["s2_own"], "deterministic")
    row(ep, "c_naive", p["c_naive"], c["c_naive_own"], "deterministic", note="Brent xatol; diff ~1e-8")
    if "c_simex_own" in c:
        mcse = 0.0005
        row(ep, "c_SIMEX (plan)", p["c_simex_plan"], c["c_simex_own"], "Monte Carlo (own seed)", tol=3 * mcse,
            note="independent SIMEX draws, B_s=50; tolerance 3 x approx MC SE")
    row(ep, "SIMEX extrapolation from project c(lambda)", p["c_simex_plan"], c.get("c_simex_from_project_c_lambda", np.nan) if "c_simex_from_project_c_lambda" in c else
        float(np.polyval(np.polyfit([float(k) for k in p["c_lambda"]], list(p["c_lambda"].values()), 2), -1)), "deterministic")
    row(ep, "D31 PB bias", p["pb_bias"], pb["pb_bias_own"], "Monte Carlo (own seed)",
        tol=3 * float(np.hypot(pb["pb_mc_se_own"], p["pb_mc_se"])), note=f"K=20 own vs K={p['pb_K']}; z={pb['z_diff']:.2f}")
    own_final = (c["c_simex_own"] if "c_simex_own" in c else p["c_simex_plan"]) - pb["pb_bias_own"]
    row(ep, "c (final = SIMEX - PB bias)", p["estimate"], own_final, "Monte Carlo",
        tol=3 * float(np.hypot(0.0005, np.hypot(pb["pb_mc_se_own"], p["pb_mc_se"]))),
        note="own SIMEX + own PB")
    row(ep, "final = c_simex_plan - pb_bias (project parts)", p["estimate"], p["c_simex_plan"] - p["pb_bias"], "arithmetic")
    # decision rule
    se = p["se_boot"]; ci = [p["estimate"] - norm.ppf(0.99) * se, p["estimate"] + norm.ppf(0.99) * se]
    dec = ("inconclusive" if (p.get("v3_false_concave") or 0) > 0.10 else "linear" if (ci[0] >= 0.9 and ci[1] <= 1.1)
           else "concave" if ci[1] < 1 else "inconclusive_convex" if ci[0] > 1 else "inconclusive")
    pc = norm.cdf((p["estimate"] - 1) / se)
    rows.append({"endpoint": ep, "quantity": "decision from CI rule", "project": p["decision"], "own": dec,
                 "absdiff": None, "kind": "rule", "pass": dec == p["decision"] and abs(ci[0] - p["ci98"][0]) < 1e-9 and abs(pc - p["p_conc"]) <= 1e-9 * max(1, pc) + 1e-18,
                 "note": f"ci98 recomputed [{ci[0]:.4f}, {ci[1]:.4f}]; p_conc {pc:.3e}"})

# E3
e3 = br["E3"]
row("E3", "k (ratio of sums)", P["E3"]["estimate"], e3["k_own"], "deterministic")
row("E3", "k internal edges only", P["E3"]["sensitivity"]["internal_edges_only"]["estimate"], e3["internal_only_own"], "deterministic")
row("E3", "sigma2_k", P["E3"]["sigma2_k"], e3["sigma2_k_own_fulltotals_ddof1"], "deterministic",
    note="v_h with HOG-total (co)variances var_sum/cov_sum_aa; diag-only b_se gives " + f"{e3['sigma2_k_own_diag_ddof1']:.3f}")
row("E3", "CATH delta log k_alpha - log k_beta", P["E3"]["cath_contrast"]["delta"], e3["cath"]["delta_own"], "deterministic")
row("E3", "n HOGs", P["E3"]["n_hogs"], br["n_hogs_own"], "exact", 0)
row("E3", "n edges", P["E3"]["n_units"], br["n_edges_own"], "exact", 0)
# E4
e4 = br["E4"]["own"]
row("E4", "beta (MixedLM REML)", P["E4"]["estimate"], e4["beta"], "deterministic (optimiser)", tol=1e-4,
    note="REML optimum; statsmodels optimiser tolerance ~1e-5")
row("E4", "mu", P["E4"]["mu"], e4["mu"], "deterministic (optimiser)", tol=1e-4)
row("E4", "gamma", P["E4"]["gamma"], e4["gamma"], "deterministic (optimiser)", tol=1e-4)
row("E4", "n cells", P["E4"]["n_units"], e4["n_cells"], "exact", 0)
row("E4", "censored fraction", P["E4"]["censored_frac"], e4["censored_frac"], "deterministic")
row("E4", "beta, plan cells (sensitivity)", P["E4"]["sensitivity"]["cells_plan_rule"]["estimate"], br["E4"]["cells_plan_rule_own"]["beta"], "deterministic (optimiser)", tol=1e-4)
# E5
row("E5", "n lineages", P["E5"]["n_lineages"], e5["n_lineages_own"], "exact", 0)
row("E5", "sigma2_obs", P["E5"]["sigma2_obs"], e5["sigma2_obs_own"], "deterministic")
row("E5", "sigma2_L = sigma2_obs - F (project F)", P["E5"]["sigma2_L"], e5["sigma2_L_own_with_project_F"], "deterministic given F",
    note="floor F is a P=100 parametric bootstrap; not re-simulated")
nfx = sum(1 for c in e5["contrast_detail"] if c["x"] is not None)
row("E5", "n contrasts (finite x_v)", P["E5"]["n_contrasts"], nfx, "exact", 0,
    note=f"{e5['n_contrasts_own']} contrasts pass the >=50-HOG rule; {nfx} have finite x_v (project's count); "
         f"{e5['n_contrasts_finite_xy_own']} have finite x_v AND y_v and enter beta")
row("E5", "beta contrast", P["E5"]["beta_contrast"], e5["beta_own_mean"], "Monte Carlo (c_v, vx_v from 200 own HOG-pigeonhole reps)",
    tol=3 * max(e5["beta_own_sd_over_seeds"], 1e-3), note=f"5 own seeds sd {e5['beta_own_sd_over_seeds']:.4f}; naive (no disattenuation) {e5['beta_naive_own']:.4f}")
# E6
k = [k for k in e6 if k.startswith("c=") and k.endswith("P1=seeded")][0]
row("E6", "E6a", P["E6"]["E6a"], e6[k]["E6a"], "deterministic (seeded P1 draw)")
row("E6", "E6b", P["E6"]["E6b"], e6[k]["E6b"], "deterministic")
row("E6", "n groups (a)", P["E6"]["n_units"]["a"], e6[k]["n_a"], "exact", 0)
row("E6", "n groups (b)", P["E6"]["n_units"]["b"], e6[k]["n_b"], "exact", 0)
kb = [k for k in e6 if k.endswith("P1=both")][0]
out = {"check": "V11", "pass_all": all(r["pass"] for r in rows), "rows": rows, "inputs": inp,
       "E6a_orientation_free": e6[kb]["E6a"],
       "notes": []}
wjson(out, OUT / "V11.json")
print(json.dumps({"pass_all": out["pass_all"], "fails": [r for r in rows if not r["pass"]]}, indent=1, default=str))
fmt = lambda v: (f"{v:.6g}" if isinstance(v, float) else str(v))
md = ["| endpoint | quantity | project | own | abs diff | tol | kind | pass | notes |", "|---|---|---|---|---|---|---|---|---|"]
for r in rows:
    md.append(f"| {r['endpoint']} | {r['quantity']} | {fmt(r['project'])} | {fmt(r['own'])} | "
              f"{'' if r['absdiff'] is None else format(r['absdiff'], '.2e')} | {r.get('tol','')} | {r['kind']} | {'PASS' if r['pass'] else 'FAIL'} | {r['note']} |")
(OUT / "V11_table.md").write_text("\n".join(md) + "\n")

# ---- Holm / E4 p re-check and notes
p4 = P["E4"]; pe4 = 2 * (1 - norm.cdf(abs(p4["estimate"]) / p4["se_boot"]))
holm = endpoint("holm")
pv = holm["p_raw"]; keys = sorted(pv, key=lambda k: pv[k]); run = 0; hadj = {}
for i, kk in enumerate(keys):
    run = max(run, min(1, (len(keys) - i) * pv[kk])); hadj[kk] = run
out["E4_wald_p_recomputed"] = pe4; out["E4_wald_p_project"] = p4["p_value"]
out["holm_recomputed"] = hadj; out["holm_match"] = all(abs(hadj[k] - holm["p_holm"][k]) <= 1e-12 + 1e-9 * hadj[k] for k in hadj)
out["notes"] = [
    "Inputs: confirm split only (5600 HOGs = split.tsv; pair tables, branch tables, E6 triplets all split=='confirm'); "
    "primary pairs 866959 = sum C(|S_h|,2) for every confirm HOG in all 3 alignments; no duplicate pair_ids; 29117 confirm "
    "E6 triplets (+7383 tune = 36500 = hogs_e6.tsv).",
    "E1/E2 n_excluded bookkeeping: the 'hog_lt_min_pairs_or_range' bucket also counts weight-0 pairs already counted in "
    "'pair_weight0_fail_or_sat' (E1: 20891 vs 20837 own; E2: 20932 vs 20877), so buckets + n_units exceed the primary-pair "
    "count by 54/55. n_units (pairs entering the fit) is correct. Reporting only.",
    "fail_frac per HOG recomputed from the pair table (not ok, or AA/metric saturated) matches hogs_endpoint_cols exactly; same 183 (tdi) / 140 (lddt) excluded HOGs.",
    "D31 wording vs implementation: DEVIATIONS says the PB pseudo-data use 'edge-level, tree-additive noise'; the code "
    "uses a tree-additive DESIGN (x_hat = X b_AA) with independent pair-level Sigma_ij noise. The independent PB here "
    "follows the implemented reading and reproduces pb_bias within MC error (E1 z=0.67, E2 z=1.06). Documentation fix suggested.",
    "2 HOGs (HOG0022006, HOG0102266) have NaN misfit_ratio and are kept in E3-E5 (NaN is not > tau). Consistent with the rule as written; worth a line in the report.",
    "E5 n_contrasts=146 counts contrasts with finite x_v; one of them has non-finite y_v (sum of 3Di path lengths <= 0 on one side) so beta uses 145. Reporting nuance.",
    "E5 floor F (P=100 parametric bootstrap needing per-HOG WLS operators) was not re-simulated; sigma2_L checked as sigma2_obs - F. beta_contrast depends on random c_v/vx_v (200 HOG-pigeonhole reps); own 5 seeds give 0.583 +- 0.008 vs project 0.572 (1.4 sd).",
    "E6a uses the plan's seeded P1 choice (stage e6_p1_choice) -> exact match; orientation-averaged E6a = %.5f (sign and conclusion unchanged)." % out["E6a_orientation_free"],
    "Branch tables (b, b_se, cov_aa, ages) were used as inputs and not re-derived from pairs; their tip correction (b = b_raw - eta/2 on terminal 3Di edges) and tau = midpoint were checked.",
]
wjson(out, OUT / "V11.json")
print("E4 p", pe4, p4["p_value"], "holm match", out["holm_match"])
