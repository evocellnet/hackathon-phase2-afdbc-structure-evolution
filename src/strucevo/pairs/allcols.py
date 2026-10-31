"""All-aligned-columns count matrices (E1/E2 sensitivity "All aligned columns", ANALYSIS_PLAN §4.1/§4.2).

The production pair tables store AA/3Di count matrices over CORE columns only (pLDDT >= 70 in both proteins;
counts_shard_*.npz keys aa/tdi); seqonly additionally stores aa_all. The all-column 3Di counts of the primary
(Foldseek) alignment are not stored, so they are recomputed here by re-running the SAME Foldseek alignment
(strucevo.pairs.align.run_foldseek, same method and flags) on the primary pairs of one split and counting AA and
3Di pairs over every aligned (M) column. Determinism check: the recomputed CORE counts must equal the stored core
counts (per pair, exact); pairs where they differ are reported and set to NaN (weight 0) in the sensitivity.

Output (cached): ROOT/data/interim/pairs_variants/all_columns/{aln}/allcols__{split}.npz
  pair_id, aa_all, tdi_all (N,20,20 int32), core_match (bool), aligned (bool)  + _meta__{split}.json
Never writes into data/interim/pairs/.
"""
from __future__ import annotations

import json
import multiprocessing as mp
import os
import shutil
import time
from pathlib import Path

import numpy as np

from . import posthoc
from .align import FS_FIELDS, cigar_to_pairs, cleanup_foldseek, run_foldseek


def _worker(task):
    posthoc._single_thread()
    st = posthoc.STORE
    out = []
    for row, qidx, tidx, qs, ts, cig in task:
        qi, ti = cigar_to_pairs(qs, ts, cig)
        _, pq, aq, sq = st.get(qidx)
        _, pt, at, stt = st.get(tidx)
        qi = np.asarray(qi, np.int64); ti = np.asarray(ti, np.int64)
        core = (pq[qi] >= posthoc.CORE_PLDDT) & (pt[ti] >= posthoc.CORE_PLDDT)
        eaq = posthoc.LUT[aq[qi]]; eat = posthoc.LUT[at[ti]]
        esq = posthoc.LUT[sq[qi]]; est = posthoc.LUT[stt[ti]]
        allm = np.ones(len(qi), bool)
        out.append((row, posthoc._counts(eaq, eat, core), posthoc._counts(esq, est, core),
                    posthoc._counts(eaq, eat, allm), posthoc._counts(esq, est, allm)))
    return out


def all_column_counts(root, aln: str, split: str, threads: int = 48, chunk: int = 60000, log=print):
    """Compute (or load cached) all-column AA/3Di counts for the PRIMARY pairs of `split` under alignment `aln`."""
    import polars as pl
    from strucevo.stats.common import ConfirmSplitGuard, confirm_allowed, guard_split
    from .store import Store
    if split == "confirm" and not confirm_allowed():
        raise ConfirmSplitGuard("all_column_counts(split='confirm') requires CONFIRM_RUN=1")
    if aln not in ("tdiaa", "tmalign"):
        raise ValueError("all-column counts are recomputed for Foldseek alignments only (seqonly stores aa_all)")
    root = Path(root)
    src = root / "data" / "interim" / "pairs" / aln
    out = root / "data" / "interim" / "pairs_variants" / "all_columns" / aln
    out.mkdir(parents=True, exist_ok=True)
    f_npz = out / f"allcols__{split}.npz"
    f_meta = out / f"_meta__{split}.json"
    sig = [[f.name, f.stat().st_size, int(f.stat().st_mtime)] for f in sorted(src.glob("shard_[0-9][0-9][0-9][0-9].parquet"))]
    if f_npz.exists() and f_meta.exists():
        meta = json.loads(f_meta.read_text())
        if meta.get("done") and meta.get("source_sig") == sig:
            return f_npz, meta
    t0 = time.time()
    # primary pairs of the split + their stored core counts (determinism check)
    ids, core_aa, core_tdi = [], [], []
    for pq in sorted(src.glob("shard_[0-9][0-9][0-9][0-9].parquet")):
        k = pq.stem[6:10]
        df = (pl.scan_parquet(pq).with_row_index("_r")
              .filter((pl.col("split") == split) & (pl.col("pair_set") == "primary")).select("_r", "pair_id", "split").collect())
        if df.height == 0:
            continue
        guard_split(df["split"].to_numpy())
        z = np.load(src / f"counts_shard_{k}.npz")
        rows = df["_r"].to_numpy()
        assert np.all(z["pair_id"][rows] == df["pair_id"].to_numpy()), f"pair_id mismatch shard {k}"
        ids.append(df["pair_id"].to_numpy()); core_aa.append(z["aa"][rows]); core_tdi.append(z["tdi"][rows])
    ids = np.concatenate(ids).astype(str); core_aa = np.concatenate(core_aa); core_tdi = np.concatenate(core_tdi)
    n = len(ids)
    store = Store(str(root / "data" / "interim" / "structures" / "store"))
    ga = [p.split("__")[0] for p in ids]; gb = [p.split("__")[1] for p in ids]
    qidx = np.array([store.idx_of_acc.get(g, -1) for g in ga]); tidx = np.array([store.idx_of_acc.get(g, -1) for g in gb])
    has = (qidx >= 0) & (tidx >= 0)
    A_aa = np.zeros((n, 20, 20), np.int32); A_t = np.zeros((n, 20, 20), np.int32)
    C_aa = np.zeros((n, 20, 20), np.int32); C_t = np.zeros((n, 20, 20), np.int32)
    aligned = np.zeros(n, bool)
    posthoc.STORE = store
    wd_root = out / f"work__{split}"
    hv_all = np.flatnonzero(has)
    with mp.get_context("fork").Pool(threads) as pool:
        for c0 in range(0, len(hv_all), chunk):
            hv = hv_all[c0:c0 + chunk]
            wd = str(wd_root / f"c{c0:08d}")
            tsv, _ = run_foldseek(store.db, store.key[qidx[hv]], store.key[tidx[hv]], wd, aln, threads=threads)
            fs = pl.read_csv(tsv, separator="\t", has_header=False, new_columns=FS_FIELDS, quote_char=None,
                             schema_overrides={"lddtfull": pl.Utf8, "cigar": pl.Utf8, "query": pl.Utf8, "target": pl.Utf8})
            row_of = {(int(qidx[i]), int(tidx[i])): int(i) for i in hv}
            fq = [store.idx_of_name[x] for x in fs["query"].to_list()]
            ft = [store.idx_of_name[x] for x in fs["target"].to_list()]
            tasks = [(row_of[(a, b)], a, b, int(qs), int(ts), cg) for a, b, qs, ts, cg in
                     zip(fq, ft, fs["qstart"].to_list(), fs["tstart"].to_list(), fs["cigar"].to_list())]
            for part in pool.imap_unordered(_worker, [tasks[i:i + 200] for i in range(0, len(tasks), 200)]):
                for row, ca, ct, aa, at in part:
                    C_aa[row] = ca; C_t[row] = ct; A_aa[row] = aa; A_t[row] = at; aligned[row] = True
            cleanup_foldseek(wd); shutil.rmtree(wd, ignore_errors=True)
            log(f"  all-columns {aln}/{split}: {min(c0 + chunk, len(hv_all))}/{len(hv_all)} pairs ({time.time() - t0:.0f}s)")
    shutil.rmtree(wd_root, ignore_errors=True)
    core_match = aligned & np.all(C_aa == core_aa, axis=(1, 2)) & np.all(C_t == core_tdi, axis=(1, 2))
    stored_nonempty = (core_aa.sum((1, 2)) > 0)
    meta = {"aln": aln, "split": split, "n_pairs": int(n), "n_aligned": int(aligned.sum()),
            "n_core_match": int(core_match.sum()),
            "n_core_mismatch_among_stored_nonempty": int((stored_nonempty & ~core_match).sum()),
            "runtime_s": round(time.time() - t0, 1), "threads": threads, "done": True, "source_sig": sig,
            "created": time.strftime("%Y-%m-%dT%H:%M:%S"),
            "note": "Foldseek alignment re-run with the production method/flags; counts over every aligned M column"}
    tmp = str(f_npz) + ".tmp.npz"
    np.savez_compressed(tmp, pair_id=ids, aa_all=A_aa, tdi_all=A_t, core_match=core_match, aligned=aligned)
    os.replace(tmp, f_npz)
    f_meta.write_text(json.dumps(meta, indent=1))
    log(f"  all-columns {aln}/{split}: done {meta}")
    return f_npz, meta
