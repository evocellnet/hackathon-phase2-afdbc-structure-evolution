"""In-house CA lDDT that reproduces Foldseek's per-residue ``lddtfull`` (ANALYSIS_PLAN §2.1, check V9).

Foldseek convention (verified in reviews/round1_algorithms_compute.md C3 and by scripts/verify/v9_lddt.py):
  * CA atoms only; the *query* is the reference structure; inclusion radius R0 = 15 Å;
  * neighbours of an aligned query residue i are *all* query residues j != i with d_q(i, j) < R0;
    neighbours that are not aligned to the target count as "not preserved";
  * a neighbour distance is preserved at threshold t if |d_q(i, j) - d_t(a(i), a(j))| < t, t in {0.5, 1, 2, 4};
  * per-residue score = mean over thresholds of (preserved / n_neighbours); residues without neighbours score 0.
``lddtfull`` lists these per-residue scores for the aligned query residues, so Foldseek's ``lddt`` field is the
mean over aligned residues (local) while the project's global LDDT is sum(per-residue) / len_query.
"""
from __future__ import annotations

import numpy as np

THRESHOLDS = (0.5, 1.0, 2.0, 4.0)


def _cdist(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    d2 = (a * a).sum(1)[:, None] + (b * b).sum(1)[None, :] - 2.0 * (a @ b.T)
    np.maximum(d2, 0.0, out=d2)
    return np.sqrt(d2)


def per_residue_lddt(xq: np.ndarray, xt: np.ndarray, qi: np.ndarray, ti: np.ndarray, r0: float = 15.0,
                     thresholds=THRESHOLDS, neighbour_mask: np.ndarray | None = None) -> np.ndarray:
    """Per-residue lDDT of the aligned query residues ``qi`` (aligned to target residues ``ti``).

    xq (Lq,3), xt (Lt,3): CA coordinates (float). qi, ti: equal-length int arrays of aligned residue indices.
    neighbour_mask: optional bool (Lq,) restricting which query residues may be neighbours.
    Returns float64 array of len(qi).
    """
    qi = np.asarray(qi, dtype=np.int64)
    ti = np.asarray(ti, dtype=np.int64)
    n = len(qi)
    if n == 0:
        return np.zeros(0)
    xq = np.asarray(xq, dtype=np.float64)
    xt = np.asarray(xt, dtype=np.float64)
    Lq = len(xq)
    dq = _cdist(xq[qi], xq)                              # (n, Lq)
    inc = dq < r0
    inc[np.arange(n), qi] = False
    if neighbour_mask is not None:
        inc &= np.asarray(neighbour_mask, bool)[None, :]
    tmap = np.full(Lq, -1, dtype=np.int64)
    tmap[qi] = ti
    aligned = tmap >= 0
    dt = _cdist(xt[ti], xt[np.where(aligned, tmap, 0)])  # (n, Lq), garbage where not aligned
    diff = np.abs(dq - dt)
    sc = np.zeros_like(diff)
    for t in thresholds:
        sc += diff < t
    sc /= len(thresholds)
    sc *= (inc & aligned[None, :])
    nn = inc.sum(1)
    out = np.where(nn > 0, sc.sum(1) / np.maximum(nn, 1), 0.0)
    return out


def lddt_scores(xq, xt, qi, ti, len_q: int | None = None, core: np.ndarray | None = None, r0: float = 15.0):
    """Convenience: (lddt_aln, lddt_global, lddt_core) from per-residue scores.

    lddt_aln = mean over aligned query residues (= Foldseek ``lddt``), lddt_global = sum / len_q,
    lddt_core = mean over aligned columns where ``core`` (bool per aligned column) is True.
    """
    pr = per_residue_lddt(xq, xt, qi, ti, r0=r0)
    len_q = len(xq) if len_q is None else len_q
    la = float(pr.mean()) if len(pr) else np.nan
    lg = float(pr.sum() / len_q)
    lc = float(pr[core].mean()) if core is not None and np.any(core) else np.nan
    return la, lg, lc
