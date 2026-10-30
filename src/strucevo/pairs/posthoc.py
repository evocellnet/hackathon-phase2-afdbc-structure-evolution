"""Per-pair post-hoc metrics from one parsed alignment (ANALYSIS_PLAN §1.2, §2.1; review C3).

For each pair (query = gene_a = shorter protein):
  core mask   aligned columns with pLDDT >= 70 in both proteins
  n_aln, n_core, cov_short = n_aln / len_a
  AA and 3Di 20x20 count matrices over core columns (AA_ORDER for both alphabets: trap #4);
  for e6_* pairs also the split halves A/B: alternating 20-residue blocks of the *query* residue index
  (block = (qi // 20) % 2; 0 -> A, 1 -> B)
  AA counts over *all* aligned columns (seqonly only; give aa_ml_seq / aa_pid_seq)
  lddt_global = sum(per-residue lDDT of aligned query residues) / len_a
  lddt_core   = mean per-residue lDDT over aligned core columns
  rmsd_core   = strucevo.distances.kabsch_rmsd over aligned core CA pairs
Per-residue lDDT comes from Foldseek ``lddtfull`` for Foldseek alignments and from ``strucevo.lddt`` for seqonly.
"""
from __future__ import annotations

import numpy as np

from ..distances import AA_ORDER, _lut, kabsch_rmsd
from ..lddt import per_residue_lddt
from .align import cigar_to_pairs, seqonly_pairs

CORE_PLDDT = 70.0
HALF_BLOCK = 20
LUT = _lut(AA_ORDER).astype(np.int64)   # same LUT for AA and 3Di letters (Q.3Di.AF is in PAML ARND order)

STORE = None  # set in the parent before forking (shared, memory-mapped)


def _counts(ca, cb, sel):
    a = ca[sel]; b = cb[sel]
    ok = (a >= 0) & (b >= 0)
    return np.bincount(a[ok] * 20 + b[ok], minlength=400).reshape(20, 20)


def pair_metrics(qidx: int, tidx: int, qi, ti, lddt_res=None, halves: bool = False, all_cols: bool = False,
                 inhouse_lddt: bool = False, store=None) -> dict:
    """Compute all post-hoc quantities for one aligned pair. ``lddt_res``: per-M-column lDDT (Foldseek) or None."""
    st = store if store is not None else STORE
    xq, pq, aq, sq = st.get(qidx)
    xt, pt, at, stt = st.get(tidx)
    qi = np.asarray(qi, np.int64); ti = np.asarray(ti, np.int64)
    len_a = len(pq)
    out = {"n_aln": len(qi)}
    core = (pq[qi] >= CORE_PLDDT) & (pt[ti] >= CORE_PLDDT)
    out["n_core"] = int(core.sum())
    out["cov_short"] = len(qi) / len_a
    eaq = LUT[aq[qi]]; eat = LUT[at[ti]]
    esq = LUT[sq[qi]]; est = LUT[stt[ti]]
    out["c_aa"] = _counts(eaq, eat, core)
    out["c_tdi"] = _counts(esq, est, core)
    if halves:
        blk = (qi // HALF_BLOCK) % 2
        for h, v in (("hA", 0), ("hB", 1)):
            m = core & (blk == v)
            out[f"c_aa_{h}"] = _counts(eaq, eat, m)
            out[f"c_tdi_{h}"] = _counts(esq, est, m)
    if all_cols:
        out["c_aa_all"] = _counts(eaq, eat, np.ones(len(qi), bool))
    if lddt_res is None:
        pr = per_residue_lddt(xq, xt, qi, ti)
    else:
        pr = np.asarray(lddt_res, float)
        if inhouse_lddt:  # V9: in-house lDDT on the same Foldseek alignment
            out["lddt_global_inhouse"] = float(per_residue_lddt(xq, xt, qi, ti).sum() / len_a)
    out["lddt_global"] = float(pr.sum() / len_a)
    out["lddt_core"] = float(pr[core].mean()) if core.any() else np.nan
    if core.sum() >= 3:
        out["rmsd_core"] = float(kabsch_rmsd(np.asarray(xq[qi[core]], np.float64), np.asarray(xt[ti[core]], np.float64)))
    else:
        out["rmsd_core"] = np.nan
    return out


def _empty(halves, all_cols):
    z = np.zeros((20, 20), np.int64)
    out = {"n_aln": 0, "n_core": 0, "cov_short": 0.0, "c_aa": z, "c_tdi": z, "lddt_global": np.nan,
           "lddt_core": np.nan, "rmsd_core": np.nan}
    if halves:
        out.update({"c_aa_hA": z, "c_aa_hB": z, "c_tdi_hA": z, "c_tdi_hB": z})
    if all_cols:
        out["c_aa_all"] = z
    return out


def _single_thread():
    """Workers run one BLAS/torch thread each (48 forked workers x 48 BLAS threads oversubscribes the node)."""
    global _LIMITED
    if not _LIMITED:
        import torch
        from threadpoolctl import threadpool_limits
        torch.set_num_threads(1)
        threadpool_limits(1)
        _LIMITED = True


_LIMITED = False


def worker_foldseek(task):
    """task: list of (row_i, qidx, tidx, qstart, tstart, cigar, lddtfull_str, halves, inhouse)."""
    _single_thread()
    res = []
    for row_i, qidx, tidx, qs, ts, cig, lf, halves, inhouse in task:
        qi, ti = cigar_to_pairs(qs, ts, cig)
        lr = np.array(lf.split(","), dtype=float) if lf else np.zeros(0)
        if len(lr) != len(qi):
            raise ValueError(f"lddtfull length {len(lr)} != M columns {len(qi)} (row {row_i})")
        res.append((row_i, pair_metrics(qidx, tidx, qi, ti, lr, halves=halves, inhouse_lddt=inhouse)))
    return res


def worker_seqonly(task):
    """task: list of (row_i, qidx, tidx, halves). Aligns with Biopython then computes metrics."""
    _single_thread()
    res = []
    for row_i, qidx, tidx, halves in task:
        qi, ti, _ = seqonly_pairs(STORE.seq(qidx), STORE.seq(tidx))
        m = pair_metrics(qidx, tidx, qi, ti, None, halves=halves, all_cols=True)
        res.append((row_i, m))
    return res
