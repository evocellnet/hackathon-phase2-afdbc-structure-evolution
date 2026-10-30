"""Unit tests for strucevo.stats (fast; < 2 min). Verification-scale checks live in scripts/verify/v*.py."""
import os
import warnings

import numpy as np
import pytest

warnings.simplefilter("ignore")

from strucevo.stats.common import (holm, guard_split, ConfirmSplitGuard, load_pairs, endpoint_record, ENDPOINT_KEYS,
                                   stage_rng)
from strucevo.stats.trees import PanelTree, prune, design_matrix, nj_splits
from strucevo.stats.branch import wls_two_pass, fit_branches
from strucevo.stats.pairwise import CurveData, CurveEngine, two_pass, fit_batch, reference_fit
from strucevo.stats.bootstrap import multinomial_weights, replicate_groups, edge_weights
from strucevo.stats.endpoints_branch import build_cells, EdgeView, fit_mixed
from strucevo.stats.reml import fast_mixed
from strucevo.stats.endpoints_misc import split_halves, interpret
from strucevo.stats.endpoints_curve import curve_decision
from strucevo.stats.synth import random_panel, make_world

NEWICK = "((((1:1,2:1)A:1,3:2)B:1,(4:1,5:1)C:2)D:1,(6:2,(7:1,8:1)E:1)F:2)R;"


def panel():
    ages = {"A": (10, 8, 12), "B": (20, 15, 25), "C": (15, 12, 18), "D": (40, 30, 50), "E": (12, 10, 14),
            "F": (30, 25, 35), "R": (100, 80, 120)}
    return PanelTree.from_newick(NEWICK, ages=ages, is_path=False)


# ------------------------------------------------------------------ guard / common
def test_confirm_guard(monkeypatch):
    monkeypatch.delenv("CONFIRM_RUN", raising=False)
    with pytest.raises(ConfirmSplitGuard):
        load_pairs("tmalign", "confirm")
    with pytest.raises(ConfirmSplitGuard):
        guard_split(np.array(["tune", "confirm"]))
    guard_split(np.array(["tune"]))
    monkeypatch.setenv("CONFIRM_RUN", "1")
    guard_split(np.array(["confirm"]))


def test_confirm_guard_in_engines(monkeypatch):
    monkeypatch.delenv("CONFIRM_RUN", raising=False)
    p = panel()
    tab = {"hog": np.array(["h"] * 3), "split": np.array(["confirm"] * 3), "species_a": np.array([1, 1, 2]),
           "species_b": np.array([2, 3, 3]), "status": np.array(["ok"] * 3), "aa_ml": np.ones(3),
           "aa_ml_se": np.ones(3), "aa_sat": np.zeros(3, bool)}
    with pytest.raises(ConfirmSplitGuard):
        fit_branches(tab, p, metrics=["aa_ml"])


def test_holm_known():
    ph = holm({"a": 0.01, "b": 0.04, "c": 0.03, "d": 0.005, "e": 0.5})
    assert ph["d"] == pytest.approx(0.025) and ph["a"] == pytest.approx(0.04)
    assert ph["c"] == pytest.approx(0.09) and ph["b"] == pytest.approx(0.09) and ph["e"] == pytest.approx(0.5)


def test_endpoint_keys():
    rec = endpoint_record(endpoint="E1")
    assert list(rec)[:len(ENDPOINT_KEYS)] == ENDPOINT_KEYS


def test_seeds_deterministic():
    assert stage_rng("bootstrap_E1", 3).random() == stage_rng("bootstrap_E1", 3).random()
    assert stage_rng("bootstrap_E1").random() != stage_rng("bootstrap_E2").random()


# ------------------------------------------------------------------ trees
def test_prune_unary_and_merge():
    p = panel()
    ht = prune(p, [1, 2, 4])          # induced: ((1,2),4) -> root D bifurcating (A-side, 4)
    assert ht.root_merged
    labels = sorted(p.label[c] for c in ht.edge_child)
    assert labels == sorted(["1", "2", "4", "A"])
    i = list(p.label[c] for c in ht.edge_child).index("A")
    assert ht.n_panel_edges[i] == 2        # A -> B -> D
    ht2 = prune(p, [1, 2, 3, 6, 7])       # root R with children D-side (B) and F-side
    assert ht2.n_params == len(ht2.edge_child) - 1


def test_design_exact_recovery():
    p = panel()
    S = [1, 2, 3, 4, 5, 6, 7, 8]
    ht = prune(p, S)
    rng = np.random.default_rng(0)
    b = rng.uniform(-0.05, 0.5, ht.n_params)       # negative lengths allowed (no clipping)
    ia, ib = np.triu_indices(len(S), 1)
    X = design_matrix(ht, ia, ib)
    d = X @ b
    f = wls_two_pass(X, d[None], np.full((1, len(d)), 0.05), np.ones((1, len(d)), bool), min_extra=0)
    assert np.allclose(f["b"][0], b, atol=1e-10)


def test_nj_matches_dendropy():
    import dendropy
    rng = np.random.default_rng(1)
    p = random_panel(12, rng=rng)
    S = p.taxid[p.is_leaf]
    ht = prune(p, S)
    ia, ib = np.triu_indices(len(S), 1)
    X = design_matrix(ht, ia, ib)
    bl = rng.uniform(0.05, 0.4, ht.n_params)
    D = np.zeros((len(S), len(S))); D[ia, ib] = X @ bl; D[ib, ia] = D[ia, ib]
    M = nj_splits(D)
    n = len(S)
    full = frozenset(range(n))

    def canon(s):
        s = frozenset(int(v) for v in s)
        return s if (len(s) < n / 2 or (len(s) == n / 2 and 0 in s)) else full - s
    ours = {canon(np.flatnonzero(M[:, e])) for e in range(M.shape[1]) if 1 < M[:, e].sum() < n - 1}
    tns = dendropy.TaxonNamespace([str(i) for i in range(len(S))])
    csv = "," + ",".join(str(i) for i in range(len(S))) + "\n" + "\n".join(
        str(i) + "," + ",".join(f"{D[i, j]:.10f}" for j in range(len(S))) for i in range(len(S)))
    pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(src=__import__("io").StringIO(csv), taxon_namespace=tns)
    t = pdm.nj_tree()
    t.encode_bipartitions()
    theirs = set()
    for nd in t.postorder_node_iter():
        if nd.is_leaf() or nd.parent_node is None:
            continue
        leaves = frozenset(int(l.taxon.label) for l in nd.leaf_iter())
        if 1 < len(leaves) < len(S) - 1:
            theirs.add(canon(leaves))
    assert ours == theirs


# ------------------------------------------------------------------ branch engine
def test_branch_tip_correction_and_hidden_paralogy():
    p = random_panel(40, rng=2)
    w = make_world(p, 5, rng=3)
    rows = {k: [] for k in ("hog", "split", "species_a", "species_b", "status", "aa_ml", "aa_ml_se", "aa_sat",
                            "tdi_ml", "tdi_ml_se", "tdi_sat")}
    for h in w.hogs:
        ia, ib, Da, D3, *_ = w.pairs_true[h]
        S = w.species[h]
        n = len(ia)
        d3 = D3 + 0.04
        Da = Da + 0.01 * np.random.default_rng(5).standard_normal(n)
        Da[n // 2] += 1.0                                 # injected hidden-paralogy outlier
        for k, v in (("hog", [h] * n), ("split", ["tune"] * n), ("species_a", S[ia]), ("species_b", S[ib]),
                     ("status", ["ok"] * n), ("aa_ml", Da), ("aa_ml_se", np.full(n, 0.01)), ("aa_sat", np.zeros(n, bool)),
                     ("tdi_ml", d3), ("tdi_ml_se", np.full(n, 0.01)), ("tdi_sat", np.zeros(n, bool))):
            rows[k].append(np.asarray(v))
    tab = {k: np.concatenate(v) for k, v in rows.items()}
    bf = fit_branches(tab, p, metrics=["aa_ml", "tdi_ml"], nu={"tdi_ml": 0.02})
    for hf in bf.hogs:
        k = len(hf.ia) // 2
        if abs(hf.fit["r"][0][k]) > 0.1:          # outlier not absorbed by a leverage-1 parameter
            assert hf.hidden_para[k]
        assert hf.hidden_para.sum() <= 4
    # 3Di: noiseless additive + eta=0.04 -> after nu=0.02 tip correction the tip-corrected sums equal truth
    E = bf.edges
    for hf in bf.hogs:
        *_, be_aa, be_3, ht = w.pairs_true[hf.hog]
        m = (E["hog"] == hf.hog) & (E["metric"] == "tdi_ml")
        assert np.isclose(E["b"][m].sum(), be_3.sum(), rtol=1e-6)


# ------------------------------------------------------------------ pairwise engine
def _curve_data(c=0.8, H=50, n=60, seed=0, noise=0.01):
    rng = np.random.default_rng(seed)
    g = np.repeat(np.arange(H), n)
    a = rng.lognormal(np.log(0.3), 0.3, H)
    x = rng.uniform(0.05, 2, len(g))
    y = 0.04 + a[g] * x ** c + noise * rng.standard_normal(len(g))
    return CurveData(g=g, x=x, y=y, sx=np.full(len(g), noise), sy=np.full(len(g), noise), w0=np.ones(len(g)),
                     eta=0.04, H=H)


def test_curve_recovers_c_and_matches_reference():
    d = _curve_data(c=0.8, noise=1e-4)
    e = CurveEngine(d)
    f = two_pass(e)
    assert abs(f["c"] - 0.8) < 2e-3
    c_ref = reference_fit(d.g, d.x, d.y, f["w"].cpu().numpy(), d.eta)
    assert abs(c_ref - f["c"]) < 1e-6
    cb, _ = fit_batch(e, w=f["w"][None].expand(3, -1))
    assert np.all(np.abs(cb - f["c"]) < 1e-6)
    cw, _ = fit_batch(e, w=f["w"][None], window=(f["c"] - 0.05, f["c"] + 0.05))
    assert abs(cw[0] - f["c"]) < 1e-6


def test_curve_decision_rules():
    assert curve_decision(1.0, 0.01)[0] == "linear"
    assert curve_decision(0.8, 0.02)[0] == "concave"
    assert curve_decision(1.3, 0.02)[0] == "inconclusive_convex"
    assert curve_decision(1.0, 0.2)[0] == "inconclusive"
    assert curve_decision(1.0, 0.01, v3_false_concave=0.2)[0] == "inconclusive"


# ------------------------------------------------------------------ bootstrap
def test_pigeonhole_weights():
    rng = np.random.default_rng(0)
    W = multinomial_weights(rng, 10)
    assert W.sum() == 10
    rows, gid = replicate_groups(np.array([2, 0, 1]), np.array([0, 0, 1, 2, 2]))
    assert len(rows) == 2 * 2 + 0 + 2 and len(np.unique(gid)) == 3
    ew = edge_weights(np.array([2.0]), np.array([1.0, 3.0, 0.0]), np.array([0, 0]), np.array([0, 2, 3]), np.array([0, 1, 2]))
    assert np.allclose(ew, [2 * 2.0, 0.0])


# ------------------------------------------------------------------ E4 cells / mixed model
def test_cell_merge_rule():
    p = panel()
    # one HOG with edges in bins 0 (tau<30), 1 (30-100): make bin-0 AA sum tiny -> merged into the next-older bin
    child = np.array([p.index["1"], p.index["2"], p.index["D"]])
    parent = np.array([p.index["A"], p.index["A"], p.index["R"]])
    ev = EdgeView(hog_code=np.zeros(3, int), hog_ids=np.array(["h"]), child=child, parent=parent,
                  is_terminal=np.array([True, True, False]), simple=np.ones(3, bool), b_aa=np.array([0.01, 0.01, 0.3]),
                  b_t=np.array([-0.01, 0.0, 0.1]), se_aa=np.ones(3), se_t=np.ones(3), cov=np.zeros(3),
                  span_ptr=np.arange(4), span_edges=child, hog_var_aa=np.ones(1), hog_var_t=np.ones(1), hog_cov=np.zeros(1),
                  n_excluded={})
    c = build_cells(ev, p.age)
    assert len(c["hog"]) == 1 and np.isclose(c["S_aa"][0], 0.32) and np.isclose(c["S_t"][0], 0.09)
    ev.b_t = np.array([-0.01, 0.0, -0.1])
    c = build_cells(ev, p.age)
    assert c["censored"][0] and c["S_t"][0] == 0.005


def test_fast_reml_matches_statsmodels():
    rng = np.random.default_rng(0)
    H = 400
    n = rng.integers(1, 6, H)
    g = np.repeat(np.arange(H), n)
    x = rng.normal(size=len(g)); z = rng.normal(size=len(g))
    u = rng.multivariate_normal([0, 0], [[0.3, 0.05], [0.05, 0.1]], H)
    y = 1 + 0.2 * x - 0.1 * z + u[g, 0] + u[g, 1] * x + rng.normal(0, 0.5, len(g))
    f = fit_mixed(y, x, z, g)
    r = fast_mixed(y, x, z, g)
    assert not f["flags"]
    assert abs(f["beta"] - r["beta"]) < 1e-5 and abs(f["gamma"] - r["gamma"]) < 1e-5


# ------------------------------------------------------------------ E6 / interpretation
def test_split_halves():
    idx = np.r_[5:30, 40:95, 100:131]
    A, B = split_halves(idx)
    assert not np.any(A & B) and np.all(A | B)
    assert np.all(A[(idx - 5) // 20 % 2 == 0])


def test_interpret_rows():
    d = {"E1": "linear", "E2": "linear", "E4": "beta_k_approx_0", "E5": "N0_accepted", "E6": "detected", "power": 0.9}
    out = interpret(d)
    assert out["row_ids"] == ["R1_linear_E4zero_N0"] and out["modifiers"] == ["duplication"]
    d["power"] = 0.5
    assert interpret(d)["row_ids"] == ["none"]
    assert "R4_E4neg" in interpret({"E1": "concave", "E4": "beta_k_lt_0", "E5": "inconclusive"})["row_ids"]
    assert "R5_N0_rejected_robust" in interpret({"E1": "linear", "E4": "beta_k_approx_0", "E5": "N0_rejected"})["row_ids"]
