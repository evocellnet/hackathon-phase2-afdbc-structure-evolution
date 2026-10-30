"""Alignments for given pair lists (ANALYSIS_PLAN §0, §1.2, trap #10).

* ``tmalign``  Foldseek TM-align (coordinate-only), primary.
* ``tdiaa``    Foldseek structurealign 3Di+AA (``--alignment-type 2``).
* ``seqonly``  Biopython PairwiseAligner: global, BLOSUM62, open -12 (BLAST 11/1 convention, C23), extend -1, end gaps free.

Foldseek runs on an exhaustive "fake prefilter" built with ``tsv2db`` from the DB keys of the given pairs,
with ``-e inf`` (structurealign) and ``--tmscore-threshold 0``; the caller asserts rows out == rows in and
labels any missing pair ``fail_align``.

Convention (verified, scripts/verify/v8_parsers.py): in Foldseek's cigar ``M`` consumes both, ``I`` consumes
the query only and ``D`` the target only; ``lddtfull`` lists one per-residue lDDT per ``M`` column.
"""
from __future__ import annotations

import os
import re
import subprocess
import time

import numpy as np

FS_FIELDS = ["query", "target", "qstart", "qend", "tstart", "tend", "qlen", "tlen", "alnlen", "fident", "cigar",
             "lddt", "alntmscore", "qtmscore", "ttmscore", "rmsd", "evalue", "bits", "lddtfull"]

_CIG = re.compile(r"(\d+)([MID])")


def cigar_to_pairs(qstart: int, tstart: int, cigar: str):
    """Foldseek cigar (1-based starts) -> (qi, ti) 0-based int arrays of aligned (M) residue pairs."""
    qi, ti = qstart - 1, tstart - 1
    Q, T = [], []
    for n, op in _CIG.findall(cigar):
        n = int(n)
        if op == "M":
            Q.append(np.arange(qi, qi + n)); T.append(np.arange(ti, ti + n))
            qi += n; ti += n
        elif op == "I":
            qi += n
        else:
            ti += n
    if not Q:
        return np.zeros(0, np.int64), np.zeros(0, np.int64)
    return np.concatenate(Q), np.concatenate(T)


def pairs_to_aln_strings(qseq: str, tseq: str, qstart: int, tstart: int, cigar: str):
    """Reconstruct Foldseek qaln/taln strings from the cigar (for check V8)."""
    qi, ti = qstart - 1, tstart - 1
    qa, ta = [], []
    for n, op in _CIG.findall(cigar):
        n = int(n)
        if op == "M":
            qa.append(qseq[qi:qi + n]); ta.append(tseq[ti:ti + n]); qi += n; ti += n
        elif op == "I":
            qa.append(qseq[qi:qi + n]); ta.append("-" * n); qi += n
        else:
            qa.append("-" * n); ta.append(tseq[ti:ti + n]); ti += n
    return "".join(qa), "".join(ta)


def run_foldseek(db: str, qkeys: np.ndarray, tkeys: np.ndarray, workdir: str, method: str, threads: int = 48,
                 extra_fields: list[str] | None = None, tmalign_fast: int = 1) -> tuple[str, dict]:
    """Align the given (query key, target key) pairs exhaustively; returns (convertalis TSV path, timings)."""
    os.makedirs(workdir, exist_ok=True)
    tm = {}
    order = np.lexsort((tkeys, qkeys))
    pre = os.path.join(workdir, "pairs.tsv")
    with open(pre, "w") as fh:
        for q, t in zip(qkeys[order], tkeys[order]):
            fh.write(f"{q}\t{t}\t0\t0\n")
    t0 = time.time()
    _run(["foldseek", "tsv2db", pre, os.path.join(workdir, "pref"), "--output-dbtype", "7", "-v", "1"])
    tm["tsv2db_s"] = time.time() - t0
    aln = os.path.join(workdir, "aln")
    t0 = time.time()
    if method == "tmalign":
        _run(["foldseek", "tmalign", db, db, os.path.join(workdir, "pref"), aln, "-a", "1", "--tmscore-threshold",
              "0", "--tmalign-fast", str(tmalign_fast), "--threads", str(threads), "-v", "1"])
    elif method == "tdiaa":
        _run(["foldseek", "structurealign", db, db, os.path.join(workdir, "pref"), aln, "-a", "1", "-e", "inf",
              "--tmscore-threshold", "0", "--lddt-threshold", "0", "--alignment-type", "2", "--threads",
              str(threads), "-v", "1"])
    else:
        raise ValueError(method)
    tm["align_s"] = time.time() - t0
    fields = FS_FIELDS + (extra_fields or [])
    out = os.path.join(workdir, "aln.tsv")
    t0 = time.time()
    _run(["foldseek", "convertalis", db, db, aln, out, "--format-output", ",".join(fields), "--threads",
          str(threads), "-v", "1"])
    tm["convertalis_s"] = time.time() - t0
    return out, tm


def cleanup_foldseek(workdir: str):
    for f in os.listdir(workdir):
        if f.startswith(("pref", "aln")) and not f.endswith(".tsv"):
            os.remove(os.path.join(workdir, f))


def _run(cmd):
    env = {k: v for k, v in os.environ.items() if k not in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                                                              "MKL_NUM_THREADS")}
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if r.returncode != 0:
        raise RuntimeError(f"{' '.join(cmd)}\n{r.stderr[-3000:]}")


# ------------------------------------------------------------------------------------------- seqonly
_ALIGNER = None
_VALID = set("ARNDCQEGHILKMFPSTWYVBZX")


def get_aligner():
    global _ALIGNER
    if _ALIGNER is None:
        from Bio.Align import PairwiseAligner, substitution_matrices
        al = PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                             open_gap_score=-12, extend_gap_score=-1)  # C23: BLAST 11/1 convention
        al.end_gap_score = 0  # end gaps free (both sequences, both ends)
        _ALIGNER = al
    return _ALIGNER


def clean_seq(s: str) -> str:
    return "".join(c if c in _VALID else "X" for c in s)


def seqonly_pairs(qseq: str, tseq: str):
    """Global BLOSUM62 alignment (first optimal path) -> (qi, ti) aligned residue pairs, score."""
    al = get_aligner()
    aln = next(iter(al.align(clean_seq(qseq), clean_seq(tseq))))
    bq, bt = aln.aligned
    if len(bq) == 0:
        return np.zeros(0, np.int64), np.zeros(0, np.int64), float(aln.score)
    qi = np.concatenate([np.arange(a, b) for a, b in bq])
    ti = np.concatenate([np.arange(a, b) for a, b in bt])
    return qi, ti, float(aln.score)
