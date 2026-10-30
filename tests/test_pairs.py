"""Tests for strucevo.pairs and strucevo.lddt (ANALYSIS_PLAN §1.2, §2.1, V5/V8/V9 unit level).

Integration tests use the demo store data/interim/pairs_demo/store (built by scripts/pairs/run_pairs.py) and are
skipped when it is absent.
"""
import json
import os

import numpy as np
import polars as pl
import pytest

from strucevo.distances import AA_ORDER, load_paml
from strucevo.lddt import per_residue_lddt
from strucevo.pairs import posthoc
from strucevo.pairs.align import cigar_to_pairs, pairs_to_aln_strings, seqonly_pairs
from strucevo.pairs.ml import Q3DI_MD5, Q3DI_PATH, lddt_D, ml_from_counts
from strucevo.pairs.pipeline import SCHEMA_ORDER, prepare_pairs

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
DEMO_STORE = os.path.join(ROOT, "data/interim/pairs_demo/store")
SEED = 20261030


# ------------------------------------------------------------------------------------------ helpers
class FakeStore:
    """Minimal in-memory stand-in for strucevo.pairs.store.Store."""

    def __init__(self, entries):
        self.entries = entries  # list of (acc, ca, plddt, aa, tdi)
        self.idx_of_acc = {e[0]: i for i, e in enumerate(entries)}
        self.len = np.array([len(e[2]) for e in entries])

    def get(self, i):
        acc, ca, pl_, aa, tdi = self.entries[i]
        return (ca, pl_, np.frombuffer(aa.encode(), np.uint8), np.frombuffer(tdi.encode(), np.uint8))

    def seq(self, i):
        return self.entries[i][3]


def random_chain(L, rng):
    steps = rng.normal(size=(L, 3))
    steps = 3.8 * steps / np.linalg.norm(steps, axis=1, keepdims=True)
    return np.cumsum(steps, 0)


def brute_lddt(xq, xt, qi, ti, r0=15.0):
    tmap = {int(q): int(t) for q, t in zip(qi, ti)}
    out = []
    for q in qi:
        nb = [j for j in range(len(xq)) if j != q and np.linalg.norm(xq[q] - xq[j]) < r0]
        if not nb:
            out.append(0.0)
            continue
        s = 0.0
        for j in nb:
            if j not in tmap:
                continue
            dq = np.linalg.norm(xq[q] - xq[j]); dt = np.linalg.norm(xt[tmap[int(q)]] - xt[tmap[j]])
            s += np.mean([abs(dq - dt) < t for t in (0.5, 1, 2, 4)])
        out.append(s / len(nb))
    return np.array(out)


# ------------------------------------------------------------------------------------------ parsers (V8 unit)
def test_cigar_parsing_convention():
    # Foldseek: M both, I query only, D target only; 1-based starts
    qi, ti = cigar_to_pairs(2, 3, "2M1I2M2D1M")
    assert qi.tolist() == [1, 2, 4, 5, 6]
    assert ti.tolist() == [2, 3, 4, 5, 8]
    qa, ta = pairs_to_aln_strings("ABCDEFGH", "abcdefghij", 2, 3, "2M1I2M2D1M")
    assert qa == "BCDEF--G" and ta == "cd-efghi"
    assert len(qa) == len(ta)


def test_cigar_empty():
    qi, ti = cigar_to_pairs(1, 1, "")
    assert len(qi) == 0 and len(ti) == 0


# ------------------------------------------------------------------------------------------ lDDT (V9 unit)
def test_lddt_identity_and_rigid_motion():
    rng = np.random.default_rng(SEED)
    x = random_chain(120, rng)
    idx = np.arange(120)
    assert np.allclose(per_residue_lddt(x, x, idx, idx), 1.0)
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    y = x @ Q.T + 5.0
    assert np.allclose(per_residue_lddt(x, y, idx, idx), 1.0)


def test_lddt_matches_bruteforce_with_partial_alignment():
    rng = np.random.default_rng(SEED + 1)
    xq = random_chain(80, rng)
    xt = random_chain(95, rng)
    xt[10:70] = xq[5:65] + rng.normal(scale=1.0, size=(60, 3))
    qi = np.arange(5, 65); ti = np.arange(10, 70)
    keep = rng.random(60) > 0.2
    qi, ti = qi[keep], ti[keep]
    a = per_residue_lddt(xq, xt, qi, ti)
    b = brute_lddt(xq, xt, qi, ti)
    assert np.allclose(a, b, atol=1e-12)
    # unaligned neighbours count as not preserved -> global LDDT of a partial alignment < 1 even if perfect
    assert per_residue_lddt(xq, xq, qi, qi).mean() < 1.0


# ------------------------------------------------------------------------------------------ seqonly
def test_seqonly_identity_and_free_end_gaps():
    s = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWUAVDIGG"
    qi, ti, sc = seqonly_pairs(s, s)
    assert (qi == ti).all() and len(qi) == len(s)
    sub = s[30:120]
    qi, ti, sc2 = seqonly_pairs(sub, s)  # end gaps free: substring aligns ungapped at its offset
    assert (ti - qi == 30).all() and len(qi) == len(sub)


def test_seqonly_gap_scores():
    from strucevo.pairs.align import get_aligner
    al = get_aligner()
    a = "WWWWWWWWWWCCCCCCCCCC"
    b = "WWWWWWWWWWGGGCCCCCCCCCC"  # 3-residue internal insertion
    # BLOSUM62 W/W = 11, C/C = 9; gap of length 3 = open(-12) + 2 * extend(-1)  (BLAST 11/1 convention, DEVIATIONS C23)
    assert al.score(a, b) == pytest.approx(10 * 11 + 10 * 9 - 12 - 2)


# ------------------------------------------------------------------------------------------ post-hoc
def _fake_pair_store(rng, L=200):
    x = random_chain(L, rng)
    pl_ = np.full(L, 90.0, np.float32); pl_[:30] = 40.0
    aa = "".join(rng.choice(list(AA_ORDER), size=L))
    tdi = "".join(rng.choice(list(AA_ORDER), size=L))
    return FakeStore([("P1", x.astype(np.float32), pl_, aa, tdi)])


def test_pair_metrics_self_pair_and_halves():
    rng = np.random.default_rng(SEED + 2)
    st = _fake_pair_store(rng)
    idx = np.arange(200)
    m = posthoc.pair_metrics(0, 0, idx, idx, None, halves=True, all_cols=True, store=st)
    assert m["n_aln"] == 200 and m["n_core"] == 170 and m["cov_short"] == 1.0
    assert m["c_aa"].sum() == 170 and np.trace(m["c_aa"]) == 170 and np.trace(m["c_tdi"]) == 170
    assert m["c_aa_all"].sum() == 200
    assert (m["c_aa_hA"] + m["c_aa_hB"] == m["c_aa"]).all()
    # halves = alternating 20-residue blocks of the query index: block 1 (20-39) is B, 40-59 is A
    m2 = posthoc.pair_metrics(0, 0, np.arange(40, 60), np.arange(40, 60), None, halves=True, store=st)
    assert m2["c_aa_hB"].sum() == 0 and m2["c_aa_hA"].sum() == 20
    assert m["lddt_global"] == pytest.approx(1.0) and m["lddt_core"] == pytest.approx(1.0)
    assert m["rmsd_core"] < 1e-6


def test_pair_metrics_uses_lddtfull_when_given():
    rng = np.random.default_rng(SEED + 3)
    st = _fake_pair_store(rng)
    qi = np.arange(50, 150)
    lf = np.full(100, 0.5)
    m = posthoc.pair_metrics(0, 0, qi, qi, lf, store=st)
    assert m["lddt_global"] == pytest.approx(50 / 200)
    assert m["cov_short"] == pytest.approx(0.5)


# ------------------------------------------------------------------------------------------ ML from counts
def test_q3di_model_file_and_order():
    import hashlib
    assert hashlib.md5(open(Q3DI_PATH, "rb").read()).hexdigest() == Q3DI_MD5
    m = load_paml(Q3DI_PATH, alphabet=AA_ORDER)
    assert m.alphabet == AA_ORDER and abs(m.pi.sum() - 1) < 1e-3


def test_ml_from_counts_floor_and_empty():
    c = np.zeros((3, 20, 20), np.int64)
    c[0] = np.diag(np.full(20, 10))          # identical, 200 sites -> d = 0, se floor 1/200
    c[2] = np.diag(np.full(20, 5)); c[2, 0, 1] = 20   # some differences
    r = ml_from_counts(c, "aa", alpha=1.0, device="cpu")
    assert r["ml"][0] == 0 and r["se"][0] == pytest.approx(1 / 200)
    assert np.isnan(r["ml"][1]) and np.isnan(r["se"][1]) and not r["sat"][1]
    assert r["ml"][2] > 0 and r["se"][2] >= 1 / 120
    r3 = ml_from_counts(c, "tdi", alpha=1.0, device="cpu")
    assert r3["ml"][0] == 0


def test_lddt_D():
    s0 = {"bin_lo": [0, 200, 400], "median": [0.20, 0.18, 0.15], "q95": [0.25, 0.22, 0.20]}
    r = lddt_D(np.array([1.0, 0.6, 0.19, 0.5, 0.14]), np.array([100, 300, 500, 100, 450]), s0)
    assert r["lddt_D"][0] == pytest.approx(0.0)
    assert r["lddt_D"][1] == pytest.approx(-np.log((0.6 - 0.18) / (1 - 0.18)))
    assert r["lddt_sat"][2] and r["lddt_D"][2] > 0          # censored but defined (S0 < S <= q95)
    assert np.isnan(r["lddt_D"][4]) and r["lddt_sat"][4]    # undefined below the floor median
    assert not r["lddt_sat"][3] and r["lddt_D"][3] > 0


# ------------------------------------------------------------------------------------------ orientation
def test_prepare_pairs_orientation_and_ids():
    st = FakeStore([("B", None, np.zeros(100), "", ""), ("A", None, np.zeros(300), "", ""),
                    ("C", None, np.zeros(100), "", "")])
    P = pl.DataFrame({"gene_a": ["A", "C", "X"], "gene_b": ["B", "B", "A"], "species_a": ["1", "2", "3"],
                      "species_b": ["4", "5", "6"], "hog": ["H1", "H1", "NA"]})
    out = prepare_pairs(P, st, 4)
    assert out["gene_a"].to_list() == ["B", "B", "X"]       # shorter first; tie B/C -> lexicographic
    assert out["gene_b"].to_list() == ["A", "C", "A"]
    assert out["species_a"].to_list() == [4, 5, 3]
    assert out["pair_id"].to_list() == ["B__A", "B__C", "X__A"]
    assert out["qidx"].to_list()[2] == -1                    # missing structure -> fail_align downstream
    assert out["shard"][0] == out["shard"][1]                # same HOG -> same shard
    with pytest.raises(ValueError):
        prepare_pairs(pl.concat([P, P]), st, 4)


# ------------------------------------------------------------------------------------------ integration
@pytest.mark.skipif(not os.path.exists(os.path.join(DEMO_STORE, "BUILD.json")), reason="demo store not built")
def test_runner_end_to_end(tmp_path):
    from strucevo.pairs.pipeline import Runner
    from strucevo.pairs.store import Store
    st = Store(DEMO_STORE)
    accs = st.meta.sort("length")["acc"].to_list()[2000:2040]
    P = pl.DataFrame({"gene_a": accs[:20] + accs[:5] + ["NOT_IN_STORE"],
                      "gene_b": accs[20:40] + accs[:5] + [accs[0]],
                      "pair_set": ["primary"] * 10 + ["e6_para"] * 10 + ["v5_self"] * 5 + ["primary"],
                      "hog": [f"H{i % 3}" for i in range(26)]})
    P = prepare_pairs(P, st, 2)
    with Runner(DEMO_STORE, str(tmp_path), threads=4, device="cpu") as R:
        infos = R.run(P)
        assert sum(i["rows_out"] for i in infos if i["aln"] == "tmalign") == 26
        for aln in ("seqonly", "tmalign", "tdiaa"):
            d = pl.read_parquet(str(tmp_path / aln / "shard_*.parquet"))
            assert d.columns == SCHEMA_ORDER and d.height == 26
            assert d.filter(pl.col("gene_a") == "NOT_IN_STORE")["status"].to_list() == ["fail_align"]
            s = d.filter(pl.col("gene_a") == pl.col("gene_b"))
            assert (s["lddt_global"] > 0.999).all() and (s["aa_ml_seq"] == 0).all()
            if aln != "seqonly":
                assert ((s["tm_short"] - 1).abs() < 1e-3).all()
            e6 = d.filter(pl.col("pair_set") == "e6_para")
            assert e6["aa_ml_hA"].is_not_nan().any() and d.filter(pl.col("pair_set") == "primary")["aa_ml_hA"].is_nan().all()
            z = np.load(sorted((tmp_path / aln).glob("counts_shard_*.npz"))[0])
            assert z["aa"].shape[1:] == (20, 20)
        # resumable: second run skips every shard
        infos2 = R.run(P)
        assert all(i["skipped"] for i in infos2)
        # a tampered shard is recomputed
        done = sorted((tmp_path / "tmalign").glob("*.done"))[0]
        info = json.load(open(done)); info["parquet_md5"] = "x"; json.dump(info, open(done, "w"))
        infos3 = R.run(P, alns=("seqonly", "tmalign"))
        assert sum(not i["skipped"] for i in infos3) == 1
