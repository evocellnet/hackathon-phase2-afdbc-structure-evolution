"""E6 (duplication positive control, ANALYSIS_PLAN §4.6), split-half helper, Holm (§5) and the §11.3
interpretation mapping (row ids only)."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.stats import norm

from .common import (Z98, endpoint_record, normal_ci, pct_ci, stage_rng, md5_arrays, holm, write_json, ROOT,
                     guard_split)
from .bootstrap import multinomial_weights

RESULTS = ROOT / "results" / "endpoints"


# ------------------------------------------------------------------------------------------ split halves
def split_halves(residue_index, block=20):
    """Alternating 20-residue blocks by the paralog's residue index: A = odd blocks (1st, 3rd, ...), B = even.
    Returns boolean masks (A, B) over the given core columns."""
    r = np.asarray(residue_index)
    blk = (r - r.min()) // block if len(r) else r
    A = (blk % 2) == 0          # block 0 is the 1st (odd-numbered) block
    return A, ~A


# ------------------------------------------------------------------------------------------ E6
def _get(t, k):
    if hasattr(t, "get_column"):
        return t.get_column(k).to_numpy()
    return np.asarray(t[k])


def run_E6(triplets, pairs, curve, B=1000, metric=("aa_ml", "tdi_ml"), input_md5=None, secondary=False):
    """triplets: columns hog, dcg_id, species, copy1, copy2, event_level, aa_pid_seq, aa_ml_seq, o1 (list), o2 (list).
    pairs: columns gene_a, gene_b, status, {x}_hA, {x}_hB, {y}_hA, {y}_hB (+ split).
    curve: {'c', 'eta', 'a_h': {hog: a_h}} from E1."""
    xm, ym = metric
    if "split" in (pairs.columns if hasattr(pairs, "columns") else pairs):
        guard_split(_get(pairs, "split"))
    ga, gb = _get(pairs, "gene_a").astype(str), _get(pairs, "gene_b").astype(str)
    st = _get(pairs, "status").astype(str)
    cols = {k: _get(pairs, k).astype(float) for k in (f"{xm}_hA", f"{xm}_hB", f"{ym}_hA", f"{ym}_hB")}
    look = {}
    for i, (a, b) in enumerate(zip(ga, gb)):
        if st[i] == "ok":
            look[(a, b)] = i; look[(b, a)] = i
    c, eta, ah = curve["c"], curve["eta"], curve["a_h"]

    def resid(g1, g2, h, xh, yh):
        i = look.get((g1, g2))
        if i is None:
            return np.nan
        x = cols[f"{xm}_h{xh}"][i]; y = cols[f"{ym}_h{yh}"][i]
        if not (np.isfinite(x) and np.isfinite(y)):
            return np.nan
        return y - (eta + ah[h] * max(x, 0.0) ** c)

    def xval(g1, g2, xh):
        i = look.get((g1, g2))
        return cols[f"{xm}_h{xh}"][i] if i is not None else np.nan

    hog = _get(triplets, "hog").astype(str); dcg = _get(triplets, "dcg_id").astype(str)
    c1 = _get(triplets, "copy1").astype(str); c2 = _get(triplets, "copy2").astype(str)
    lev = _get(triplets, "event_level").astype(str)
    pid = _get(triplets, "aa_pid_seq").astype(float); dseq = _get(triplets, "aa_ml_seq").astype(float)
    O1 = list(_get(triplets, "o1")); O2 = list(_get(triplets, "o2"))
    rng = stage_rng("e6_p1_choice", 0)
    flip = rng.random(len(hog)) < 0.5
    excl = {"leaf_level": 0, "pid_gt_0.95_or_dseq_gt_2": 0, "no_E1_fit": 0, "no_O1": 0, "no_O2": 0, "missing_pairs_a": 0,
            "missing_pairs_b": 0}
    checks = {"o1_o2_disjoint": True, "no_leaf_events_used": True}
    da, db = {}, {}
    for i in range(len(hog)):
        if lev[i] == "leaf":
            excl["leaf_level"] += 1; continue
        if not (pid[i] <= 0.95 and dseq[i] <= 2):
            excl["pid_gt_0.95_or_dseq_gt_2"] += 1; continue
        if hog[i] not in ah or not np.isfinite(ah[hog[i]]):
            excl["no_E1_fit"] += 1; continue
        o1 = list(O1[i]) if O1[i] is not None else []
        o2 = list(O2[i]) if O2[i] is not None else []
        if set(o1) & set(o2):
            checks["o1_o2_disjoint"] = False
        if len(o1) < 1:
            excl["no_O1"] += 1; continue
        p1, p2 = (c2[i], c1[i]) if flip[i] else (c1[i], c2[i])
        h = hog[i]
        rpp = resid(p1, p2, h, "A", "B")
        ro = [resid(p1, o, h, "A", "B") for o in o1]
        ro = [v for v in ro if np.isfinite(v)]
        if np.isfinite(rpp) and ro:
            da.setdefault((h, dcg[i]), []).append(rpp - np.mean(ro))
        else:
            excl["missing_pairs_a"] += 1
        if len(o2) < 1:
            excl["no_O2"] += 1; continue
        m1 = np.nanmean([xval(c1[i], o, "A") for o in o1]); m2 = np.nanmean([xval(c2[i], o, "A") for o in o1])
        if not (np.isfinite(m1) and np.isfinite(m2)):
            excl["missing_pairs_b"] += 1; continue
        mdo, ldo = (c1[i], c2[i]) if m1 > m2 else (c2[i], c1[i])
        ex = [resid(mdo, o, h, "B", "B") - resid(ldo, o, h, "B", "B") for o in o2]
        ex = [v for v in ex if np.isfinite(v)]
        if ex:
            db.setdefault((h, dcg[i]), []).append(np.mean(ex))
        else:
            excl["missing_pairs_b"] += 1

    def agg(dd):
        keys = list(dd)
        hs = np.array([k[0] for k in keys]); vals = np.array([np.mean(dd[k]) for k in keys])
        return hs, vals
    ha, va = agg(da); hb, vb = agg(db)
    uh = np.unique(np.r_[ha, hb]) if len(ha) + len(hb) else np.array([])
    pos = {h: k for k, h in enumerate(uh)}
    ia = np.array([pos[h] for h in ha], int); ib = np.array([pos[h] for h in hb], int)
    E6a = float(np.mean(va)) if len(va) else np.nan
    E6b = float(np.mean(vb)) if len(vb) else np.nan
    brng = stage_rng("bootstrap_E6", 0)
    ra, rb = [], []
    for _ in range(B):
        W = multinomial_weights(brng, len(uh))
        wa, wb = W[ia], W[ib]
        ra.append(np.sum(wa * va) / wa.sum() if wa.sum() > 0 else np.nan)
        rb.append(np.sum(wb * vb) / wb.sum() if wb.sum() > 0 else np.nan)
    sa, sb = float(np.nanstd(ra, ddof=1)), float(np.nanstd(rb, ddof=1))
    pa = float(1 - norm.cdf(E6a / sa)) if sa > 0 else np.nan
    pb = float(1 - norm.cdf(E6b / sb)) if sb > 0 else np.nan
    p = float(min(1.0, 2 * np.nanmin([pa, pb]))) if np.isfinite([pa, pb]).any() else np.nan
    rec = endpoint_record(
        endpoint="E6" + ("_secondary_lddtD" if secondary else ""),
        estimand="(a) E6a = mean_dcg[r(P1,P2) - mean_O1 r(P1,O)] (x half A, y half B); (b) E6b = mean_dcg MDO-LDO excess (half B)",
        estimate=E6a, se_boot=sa, ci98=normal_ci(E6a, sa), ci_pct98=pct_ci(ra),
        one_sided_bound={"E6a_lower99": E6a - Z98 * sa, "E6b_lower99": E6b - Z98 * sb},
        p_value=p, p_holm=None, decision="pending_holm", n_hogs=len(uh), n_units={"a": len(va), "b": len(vb)},
        n_excluded=excl, sensitivity={}, alignment_robust=None,
        flags=(["not_split_half_valid"] if secondary else []), input_md5=input_md5 or md5_arrays(va, vb),
        E6a=E6a, E6b=E6b, se_a=sa, se_b=sb, ci98_b=normal_ci(E6b, sb), ci_pct98_b=pct_ci(rb), p_a=pa, p_b=pb,
        replication_direction={"a": float(np.sign(E6a)) if np.isfinite(E6a) else None,
                               "b": float(np.sign(E6b)) if np.isfinite(E6b) else None,
                               "note": "compare with the direction reported in PMID 40813247 before interpreting"},
        checks=checks, B=B)
    return rec


# ------------------------------------------------------------------------------------------ Holm + interpretation
HOLM_SET = ("E1", "E2", "E4", "E5", "E6")


def apply_holm(results_dir=RESULTS):
    d = Path(results_dir)
    recs = {e: json.loads((d / f"{e}.json").read_text()) for e in HOLM_SET if (d / f"{e}.json").exists()}
    pv = {e: recs[e]["p_value"] if e in recs else None for e in HOLM_SET}
    ph = holm(pv)
    for e, r in recs.items():
        r["p_holm"] = ph[e]
        if e == "E4":
            if ph[e] < 0.05:
                r["decision"] = "beta_k_gt_0" if -r["estimate"] > 0 else "beta_k_lt_0"
            else:
                r["decision"] = "beta_k_approx_0"
        if e == "E6":
            r["decision"] = "detected" if ph[e] < 0.05 else "not_detected"
        write_json(r, d / f"{e}.json")
    out = {"alpha": 0.05, "m": len(HOLM_SET), "p_raw": pv, "p_holm": ph,
           "missing": [e for e in HOLM_SET if e not in recs]}
    write_json(out, d / "holm.json")
    return out


ROWS = ["R1_linear_E4zero_N0", "R2_concave_E4pos_N0", "R3_linear_E4pos_robust", "R4_E4neg",
        "R5_N0_rejected_robust", "R6_artefact_suspect"]


def interpret(decisions: dict) -> dict:
    """decisions: {'E1': class, 'E2': class, 'E4': 'beta_k_gt_0'|'beta_k_lt_0'|'beta_k_approx_0',
    'E4_robust_alpha_timescale': bool, 'E5': 'N0_accepted'|'N0_rejected'|'artefact_suspect'|'inconclusive',
    'E6': 'detected'|'not_detected', 'E6_positive': bool, 'power': float}. Returns row ids only."""
    e1, e2 = decisions.get("E1"), decisions.get("E2")
    curv = e1 if (e2 in (e1, "inconclusive", None)) else "discordant"
    e4 = decisions.get("E4"); e5 = decisions.get("E5")
    power = decisions.get("power")
    if e5 == "N0_accepted" and (power is None or power < 0.80):
        e5 = "inconclusive"
    rows = []
    n0 = e5 == "N0_accepted"
    if curv == "linear" and e4 == "beta_k_approx_0" and n0:
        rows.append(ROWS[0])
    if curv == "concave" and e4 == "beta_k_gt_0" and n0:
        rows.append(ROWS[1])
    if curv == "linear" and e4 == "beta_k_gt_0" and decisions.get("E4_robust_alpha_timescale"):
        rows.append(ROWS[2])
    if e4 == "beta_k_lt_0":
        rows.append(ROWS[3])
    if e5 == "N0_rejected":
        rows.append(ROWS[4])
    if e5 == "artefact_suspect":
        rows.append(ROWS[5])
    mods = []
    if decisions.get("E6") == "detected" and decisions.get("E6_positive", True) and n0:
        mods.append("duplication")
    return {"row_ids": rows or ["none"], "modifiers": mods}


def write_interpretation(results_dir=RESULTS):
    d = Path(results_dir)
    r = {e: json.loads((d / f"{e}.json").read_text()) for e in ("E1", "E2", "E4", "E5", "E6") if (d / f"{e}.json").exists()}
    pw = None
    if "E5" in r:
        pw = (r["E5"].get("power") or {}).get("f=1.5", {}).get("detect_rate")
    dec = {"E1": r.get("E1", {}).get("decision"), "E2": r.get("E2", {}).get("decision"),
           "E4": r.get("E4", {}).get("decision"),
           "E4_robust_alpha_timescale": r.get("E4", {}).get("robust_alpha_timescale", False),
           "E5": r.get("E5", {}).get("decision"), "E6": r.get("E6", {}).get("decision"),
           "E6_positive": (r.get("E6", {}).get("E6a") or 0) > 0 or (r.get("E6", {}).get("E6b") or 0) > 0, "power": pw}
    out = interpret(dec)
    write_json(out, d / "interpretation.json")
    return out


def e4_robust_alpha_timescale(e4_rec, names_alpha=("alpha_x1.5", "alpha_div1.5"), names_age=("ages_strassert2021", "ages_betts2018")):
    """Same sign and unadjusted p < 0.05 in both the alpha and the deep-age sensitivities (§4.4)."""
    s = e4_rec.get("sensitivity", {})
    sign = np.sign(e4_rec["estimate"])
    need = list(names_alpha) + list(names_age)
    if not all(n in s for n in need):
        return False
    return all(np.sign(s[n]["estimate"]) == sign and s[n].get("p_unadjusted", 1) < 0.05 for n in need)
