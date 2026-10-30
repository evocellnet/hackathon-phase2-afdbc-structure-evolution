"""Sequence-only selection distances aa_ml_seq / aa_pid_seq, computed with exactly the code path of the pairs
pipeline's `seqonly` stage (strucevo.pairs.align.seqonly_pairs -> AA counts over all aligned columns with
AA_ORDER -> strucevo.pairs.ml.ml_from_counts(LG+G4, d_max=20)), but without the structure store.

Orientation follows the pipeline: query = shorter protein, ties broken lexicographically by accession.
Validated against data/interim/pairs_demo seqonly tables (see tests/test_hogs.py).
"""
from __future__ import annotations

import gzip
import multiprocessing as mp
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

from strucevo.distances import AA_ORDER, _lut  # noqa: E402
from strucevo.pairs.align import seqonly_pairs  # noqa: E402

LUT = _lut(AA_ORDER).astype(np.int64)
SEQ: dict = {}           # acc -> sequence (set in the parent before forking)


def read_fasta_dir(d: Path, accs=None) -> dict:
    out = {}
    for f in sorted(Path(d).glob("*.fasta*")):
        op = gzip.open if str(f).endswith(".gz") else open
        acc, buf = None, []
        with op(f, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    if acc is not None and (accs is None or acc in accs):
                        out[acc] = "".join(buf)
                    h = line[1:].split()[0]
                    acc = h.split("|")[1] if "|" in h else h
                    buf = []
                else:
                    buf.append(line.strip())
        if acc is not None and (accs is None or acc in accs):
            out[acc] = "".join(buf)
    return out


def orient(a: str, b: str):
    la, lb = len(SEQ[a]), len(SEQ[b])
    return (b, a) if (la > lb or (la == lb and a > b)) else (a, b)


def _counts_one(pair):
    a, b = orient(*pair)
    qa, tb = SEQ[a], SEQ[b]
    qi, ti, _ = seqonly_pairs(qa, tb)
    ea = LUT[np.frombuffer(qa.encode(), np.uint8)][qi]
    eb = LUT[np.frombuffer(tb.encode(), np.uint8)][ti]
    ok = (ea >= 0) & (eb >= 0)
    return np.bincount(ea[ok] * 20 + eb[ok], minlength=400).reshape(20, 20).astype(np.int32)


def _worker(chunk):
    return [_counts_one(p) for p in chunk]


def seq_distances(pairs: list, threads: int = 8, alpha: float = 1.0, chunk: int = 200, pool=None):
    """pairs: list of (acc_a, acc_b). Returns (aa_ml_seq, aa_pid_seq, sat) arrays (NaN if no aligned columns)."""
    from strucevo.pairs.ml import ml_from_counts
    if not pairs:
        e = np.zeros(0)
        return e, e, np.zeros(0, bool)
    chunks = [pairs[i:i + chunk] for i in range(0, len(pairs), chunk)]
    if threads <= 1 and pool is None:
        res = [_worker(c) for c in chunks]
    else:
        own = pool is None
        if own:
            pool = mp.get_context("fork").Pool(threads)
        res = pool.map(_worker, chunks, chunksize=1)
        if own:
            pool.close(); pool.join()
    C = np.stack([c for r in res for c in r])
    r = ml_from_counts(C, "aa", alpha)
    return r["ml"], 1.0 - r["p"], r["sat"]
