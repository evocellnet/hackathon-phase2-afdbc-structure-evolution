"""Unit tests for strucevo.distances (fast; no external tools)."""
import json
import os

import numpy as np
import pytest
import torch
from conftest import GLOBAL_SEED, RESULTS

from strucevo.distances import (AA_ORDER, THREEDI_ORDER, RateModel, counting_gtr, encode, estimate_gtr,
                                gamma_rates, kabsch_rmsd, load_model, load_paml, ml_distance, neg_log,
                                one_minus, p_distance, pair_counts, plddt_mask, poisson_distance,
                                kimura_distance, simple_distances, simulate_counts)

DEV = "cuda" if torch.cuda.is_available() else "cpu"


# ---------------------------------------------------------------- models / parsing
@pytest.mark.parametrize("name", ["LG", "WAG"])
def test_model_properties(name):
    m = load_model(name)
    Q = m.Q()
    assert np.allclose(Q.sum(1), 0, atol=1e-12)
    assert np.isclose(-(m.pi * np.diag(Q)).sum(), 1.0)
    F = m.pi[:, None] * Q                      # detailed balance
    assert np.allclose(F, F.T, atol=1e-12)
    P = m.P(0.7)
    assert np.allclose(P.sum(1), 1)
    assert np.allclose(m.pi @ P, m.pi, atol=1e-10)


def test_lg_known_values():
    m = load_model("LG")
    assert np.isclose(m.pi[0], 0.079066) and np.isclose(m.pi[-1], 0.069147)
    assert np.isclose(m.S[1, 0], 0.425093)    # LG A<->R exchangeability in lg.dat
    assert m.alphabet == AA_ORDER


def test_paml_roundtrip(tmp_path):
    m = load_model("WAG")
    p = tmp_path / "wag_rt.dat"
    m.to_paml(str(p))
    m2 = load_paml(str(p), AA_ORDER)
    assert np.allclose(m2.Q(), m.Q(), atol=1e-5)


def test_gamma_rates():
    r, w = gamma_rates(1.0, 4)
    # known Yang (1994) mean rates for alpha=1, K=4 (also printed by AliSim/IQ-TREE)
    assert np.allclose(r, [0.1369538, 0.4767519, 1.0000000, 2.3862944], atol=1e-6)
    assert np.isclose((r * w).sum(), 1)
    r, w = gamma_rates(0.5, 4, pinv=0.2)
    assert np.isclose((r * w).sum(), 1) and np.isclose(w.sum(), 0.8)


# ---------------------------------------------------------------- counts and simple distances
def test_pair_counts_gaps_and_mask():
    a = "AC-DEX"
    b = "ACG-EE"
    C = pair_counts(a, b)
    assert C.sum() == 3                          # A/A, C/C, E/E ; gap columns and X excluded
    s = simple_distances(C)
    assert s["p"][0] == 0
    m = np.array([1, 0, 1, 1, 1, 1], bool)
    assert pair_counts(a, b, masks=m).sum() == 2
    assert p_distance("AAAA", "AAAC") == pytest.approx(0.25)
    assert poisson_distance("AAAA", "AAAC") == pytest.approx(-np.log(0.75))
    assert kimura_distance("AAAA", "AAAC") == pytest.approx(-np.log(1 - 0.25 - 0.2 * 0.0625))
    # Kimura undefined for p > ~0.854 -> NaN + flag
    s = simple_distances(pair_counts("AAAAAAAAAA", "CCCCCCCCCA"))
    assert np.isnan(s["kimura"][0]) and s["kimura_saturated"][0]


def test_plddt_mask():
    m = plddt_mask([90, 60, np.nan, 80], [75, 90, 90, 69.9])
    assert m.tolist() == [True, False, False, False]


def test_encode_3di_alphabet():
    e = encode("ACDy-X", THREEDI_ORDER)
    assert e.tolist() == [0, 1, 2, 19, -1, -1]


# ---------------------------------------------------------------- ML distance: controls
def test_identity_gives_zero():
    m = load_model("LG")
    rng = np.random.default_rng(GLOBAL_SEED)
    seqs = ["".join(rng.choice(list(AA_ORDER), 300, p=m.pi)) for _ in range(20)]
    r = ml_distance(pair_counts(seqs, seqs), m, alpha=1.0, device=DEV)
    assert np.all(r.d == 0) and np.all(r.at_zero) and not r.saturated.any()
    assert np.all(np.isnan(r.se))


def test_empty_pair_is_nan_and_flagged():
    m = load_model("LG")
    r = ml_distance(pair_counts("A-C", "-A-"), m, device=DEV)
    assert np.isnan(r.d[0]) and r.saturated[0] and r.n_sites[0] == 0


def test_random_unrelated_flagged_saturated():
    m = load_model("LG")
    rng = np.random.default_rng(GLOBAL_SEED + 1)
    out = {}
    for L in (100, 300, 1000):
        a = rng.choice(20, (1000, L), p=m.pi)
        b = rng.choice(20, (1000, L), p=m.pi)
        C = np.zeros((1000, 400), np.int64)
        np.add.at(C, (np.arange(1000)[:, None], a * 20 + b), 1)
        r = ml_distance(C.reshape(1000, 20, 20), m, alpha=1.0, device=DEV)
        out[L] = dict(frac_saturated=float(r.saturated.mean()), frac_at_cap=float(r.sat_cap.mean()),
                      frac_identity_le_random=float(r.sat_identity.mean()))
        assert r.saturated.mean() >= 0.95, (L, out[L])
    json.dump(out, open(os.path.join(RESULTS, "random_pairs.json"), "w"), indent=1)


def test_related_not_flagged():
    m = load_model("LG")
    C = simulate_counts(m, np.full(500, 3.0), 300, 1.0, rng=GLOBAL_SEED + 2)
    r = ml_distance(C, m, alpha=1.0, device=DEV)
    assert r.saturated.mean() < 0.01


def test_ml_matches_scipy_reference():
    """Independent check: maximise the same likelihood with scipy expm + bounded scalar optimiser."""
    from scipy.linalg import expm
    from scipy.optimize import minimize_scalar
    m = load_model("WAG")
    rr, ww = gamma_rates(0.7, 4, pinv=0.1)
    C = simulate_counts(m, [0.03, 0.4, 1.3, 2.5], 400, 0.7, pinv=0.1, rng=GLOBAL_SEED + 3)
    r = ml_distance(C, m, alpha=0.7, pinv=0.1, device=DEV)
    Q = m.Q()
    for n in range(len(C)):
        def nll(d):
            P = 0.1 * np.eye(20) + sum(w * expm(Q * ri * d) for ri, w in zip(rr, ww))
            return -(C[n] * np.log(P)).sum()
        ref = minimize_scalar(nll, bounds=(1e-6, 10), method="bounded", options=dict(xatol=1e-9)).x
        assert r.d[n] == pytest.approx(ref, rel=1e-4, abs=1e-6)


def test_se_matches_empirical_sd():
    m = load_model("LG")
    d = np.repeat([0.1, 0.5, 1.5], 2000)
    C = simulate_counts(m, d, 500, 1.0, rng=GLOBAL_SEED + 4)
    r = ml_distance(C, m, alpha=1.0, device=DEV)
    for x in (0.1, 0.5, 1.5):
        s = d == x
        assert np.nanmean(r.se[s]) == pytest.approx(r.d[s].std(), rel=0.08)


def test_pinv_path_unbiased():
    m = load_model("LG")
    d = np.repeat([0.2, 1.0, 2.0], 1000)
    C = simulate_counts(m, d, 1000, 0.8, pinv=0.15, rng=GLOBAL_SEED + 5)
    r = ml_distance(C, m, alpha=0.8, pinv=0.15, device=DEV)
    for x in (0.2, 1.0, 2.0):
        assert r.d[d == x].mean() / x == pytest.approx(1, abs=0.03)


# ---------------------------------------------------------------- 3Di path / GTR estimation
def _random_gtr(rng, alphabet=THREEDI_ORDER):
    K = len(alphabet)
    S = np.exp(rng.normal(0, 1.0, (K, K)))
    S = np.triu(S, 1)
    S = S + S.T
    pi = rng.dirichlet(np.full(K, 5.0))
    return RateModel(S, pi, alphabet, "synthetic-3Di", "random test model")


def test_3di_alphabet_ml_distance_from_strings():
    """General path on the 3Di alphabet: simulate sequences in Python under a synthetic GTR, recover d."""
    rng = np.random.default_rng(GLOBAL_SEED + 6)
    m = _random_gtr(rng)
    r_, w_ = gamma_rates(1.0, 4)
    letters = np.array(list(THREEDI_ORDER))
    L, R = 1000, 200
    out = {}
    for d in (0.1, 0.5, 1.0, 2.0):
        A, B = [], []
        for _ in range(R):
            a = rng.choice(20, L, p=m.pi)
            cat = rng.choice(4, L)
            b = np.empty(L, int)
            for k in range(4):
                P = m.P(d * r_[k])
                idx = np.nonzero(cat == k)[0]
                cum = P[a[idx]].cumsum(1)
                b[idx] = (rng.random(len(idx))[:, None] > cum).sum(1).clip(max=19)
            A.append("".join(letters[a])); B.append("".join(letters[b]))
        res = ml_distance(pair_counts(A, B, THREEDI_ORDER), m, alpha=1.0, device=DEV)
        out[d] = float(res.d.mean())
        assert res.d.mean() / d == pytest.approx(1, abs=0.03), (d, res.d.mean())


def test_estimate_gtr_recovers_model():
    rng = np.random.default_rng(GLOBAL_SEED + 7)
    true = _random_gtr(rng)
    d = rng.uniform(0.05, 0.6, 4000)
    C = simulate_counts(true, d, 300, 1.0, rng=rng)
    est = estimate_gtr(C, THREEDI_ORDER, alpha=1.0, device=DEV)
    cnt = counting_gtr(C, THREEDI_ORDER)
    iu = np.triu_indices(20, 1)
    ls_true = np.log(true.normalized_S()[iu])
    r_ml = np.corrcoef(ls_true, np.log(est.S[iu]))[0, 1]
    r_ct = np.corrcoef(ls_true, np.log(cnt.S[iu]))[0, 1]
    # held-out pairs: distance under estimated model vs truth
    dt = np.repeat([0.1, 0.5, 1.0, 2.0], 500)
    Ct = simulate_counts(true, dt, 500, 1.0, rng=rng)
    de = ml_distance(Ct, est, alpha=1.0, device=DEV).d
    dc = ml_distance(Ct, cnt, alpha=1.0, device=DEV).d
    res = dict(corr_logS_ml=float(r_ml), corr_logS_counting=float(r_ct),
               rel_bias_ml={x: float(de[dt == x].mean() / x - 1) for x in (0.1, 0.5, 1.0, 2.0)},
               rel_bias_counting={x: float(dc[dt == x].mean() / x - 1) for x in (0.1, 0.5, 1.0, 2.0)})
    json.dump(res, open(os.path.join(RESULTS, "estimate_gtr.json"), "w"), indent=1)
    assert r_ml > 0.9, res
    for x in (0.1, 0.5, 1.0, 2.0):
        assert abs(res["rel_bias_ml"][x]) < 0.05, res


# ---------------------------------------------------------------- transforms
def test_transforms():
    assert np.allclose(one_minus([1.0, 0.3]), [0.0, 0.7])
    assert np.allclose(neg_log([1.0, np.exp(-2)]), [0.0, 2.0])
    assert np.isfinite(neg_log(0.0))


# ---------------------------------------------------------------- Kabsch
def _rand_rot(rng, n):
    q = rng.normal(size=(n, 4))
    q /= np.linalg.norm(q, axis=1, keepdims=True)
    w, x, y, z = q.T
    return np.stack([
        np.stack([1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)], -1),
        np.stack([2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)], -1),
        np.stack([2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)], -1)], 1)


def test_kabsch_rigid_motion_zero():
    rng = np.random.default_rng(GLOBAL_SEED + 8)
    Y = rng.normal(0, 10, (64, 150, 3))
    R = _rand_rot(rng, 64)
    t = rng.normal(0, 50, (64, 1, 3))
    X = np.einsum("bij,blj->bli", R, Y) + t
    rmsd, Rf, tf = kabsch_rmsd(X, Y, return_transform=True)
    assert np.all(rmsd < 1e-8)
    assert np.allclose(np.einsum("bij,blj->bli", Rf, X) + tf[:, None], Y, atol=1e-6)
    assert np.allclose(np.linalg.det(Rf), 1)


def test_kabsch_reflection_not_zero():
    rng = np.random.default_rng(GLOBAL_SEED + 9)
    Y = rng.normal(0, 10, (100, 3))
    X = Y * np.array([-1, 1, 1])
    assert kabsch_rmsd(X, Y) > 1.0


def test_kabsch_noise_expected_rmsd_and_mask():
    rng = np.random.default_rng(GLOBAL_SEED + 10)
    B, L, sig = 200, 500, 0.5
    Y = rng.normal(0, 10, (B, L, 3))
    E = rng.normal(0, sig, (B, L, 3))
    R = _rand_rot(rng, B)
    X = np.einsum("bij,blj->bli", R, Y + E) + 3.0
    rmsd = kabsch_rmsd(X, Y)
    expected = sig * np.sqrt((3 * L - 6) / L)      # 6 rigid-body dof absorbed
    assert np.mean(rmsd) == pytest.approx(expected, rel=0.01)
    # compare against an independent per-structure numpy reference implementation
    for b in range(5):
        Hc = (Y[b] - Y[b].mean(0)).T @ (X[b] - X[b].mean(0))
        U, S, Vt = np.linalg.svd(Hc)
        dd = np.sign(np.linalg.det(Vt.T @ U.T))
        Rr = Vt.T @ np.diag([1, 1, dd]) @ U.T          # rotates Y onto X
        ref = np.sqrt((((Y[b] - Y[b].mean(0)) @ Rr.T - (X[b] - X[b].mean(0))) ** 2).sum(1).mean())
        assert rmsd[b] == pytest.approx(ref, rel=1e-9)
    # mask: corrupt 50 positions massively, mask them out -> same as RMSD on the rest
    Xc = X.copy()
    Xc[:, :50] += 100.0
    mask = np.ones((B, L), bool)
    mask[:, :50] = False
    r_mask = kabsch_rmsd(Xc, Y, mask=mask)
    r_sub = kabsch_rmsd(X[:, 50:], Y[:, 50:])
    assert np.allclose(r_mask, r_sub, atol=1e-9)
    # torch / GPU path gives the same answer
    if torch.cuda.is_available():
        rt = kabsch_rmsd(torch.tensor(X, device="cuda"), torch.tensor(Y, device="cuda")).cpu().numpy()
        assert np.allclose(rt, rmsd, atol=1e-9)


def test_kabsch_known_two_cluster_perturbation():
    """Analytic case: symmetric 6-atom cross, the +x / -x arms stretched by e. Centroid and optimal rotation are
    unchanged (identity), so RMSD = sqrt(2 e^2 / 6)."""
    Y = np.array([[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float) * 3
    e = 0.3
    X = Y.copy()
    X[0, 0] += e
    X[1, 0] -= e          # symmetric stretch: centroid unchanged, optimal rotation = identity
    assert kabsch_rmsd(X, Y) == pytest.approx(np.sqrt(2 * e * e / 6), rel=1e-9)
