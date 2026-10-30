#!/usr/bin/env python3
"""ANALYSIS_PLAN §2.3 SE calibration: block column bootstrap of the boot_sub pairs, per alignment.

Inputs (ROOT-relative): data/interim/pairs_bootsub/{aln}/shard_*.parquet + counts_shard_*.npz and the per-column
dumps data/interim/pairs_bootsub/columns/{aln}/cols_shard_*.npz written by
  run_pairs.py --pairs data/interim/pairs_bootsub/bootsub_pairs.tsv --out data/interim/pairs_bootsub \
               --dump-columns data/interim/pairs_bootsub/columns      (list from select_bootsub.py)
Tuned alpha from --alpha-aa/--alpha-tdi, else config/tuning.yaml, else results/tuning/alpha.json. S0 per C28.

Per pair (status ok), B = 100 circular block-bootstrap replicates over the aligned columns in query order:
block length 20 columns; ceil(n/20) blocks with uniform random start (wrapping around the end) are concatenated
and truncated to the original n_aln columns. For each replicate: core (pLDDT >= 70 both) AA and 3Di 20x20
counts -> ML distances at the tuned alpha (GPU, ml_distance, d_max 20); S* = sum per-residue lDDT of the
resampled columns / len_a (lddt_global analogue); D* = lddt_D(S*) with the C28 S0 table.
  kappa_m(bin)  = median over pairs of SD_boot(d_m*) / Fisher SE_m (SE at the tuned alpha, interior pairs:
                  not saturated, d > 0, <= 10% capped replicates), m in {aa_ml, tdi_ml}; bins = 10 aa_ml bins
                  (the pair's own aa_ml at tuned alpha in that alignment).
  rho(bin)      = Fisher-z mean over pairs of the per-pair correlation of the replicate deviations
                  (aa_ml*, tdi_ml*) -> aa_tdi; (aa_ml*, lddt_D*) -> aa_lddt.
  lddt SE model = median over pairs of SD_boot(S*) by length bin (floors bins) x S bin;
                  lddt_D_se = sd_S(L, S) / (S - S0(L))  (delta method, |dD/dS| = 1/(S - S0)).
Bins with < 20 eligible pairs (cells: < 5) count as empty and are filled from the nearest non-empty bin
(flagged *_filled; n per bin reported; raw values for bins with >= 1 pair are in details).

  python scripts/tuning/se_calibration.py [--root R] [--alns tdiaa seqonly tmalign] [--B 100]
Output ROOT/results/tuning/se_calibration.json. Top-level keys = config/tuning.yaml nesting:
  kappa.{aln}.{bin_edges, aa_ml, tdi_ml}; rho.{aln}.{bin_edges, aa_tdi, aa_lddt};
  lddt_se_model.{aln}.{len_bin_lo, s_bin_lo, sd_S}   (sd_S[len_bin][s_bin]); details under "details".
"""
from __future__ import annotations

import argparse
import glob
import time
from pathlib import Path

from se_common import (ALNS, LEN_BIN_LO, add_root_arg, bin_index, fill_nearest, np, out_json, pl, provenance,
                       read_alpha, read_s0, split_md5, torch, write_json)

from strucevo.distances import ml_distance
from strucevo.pairs.ml import D_MAX, lddt_D, ml_from_counts, models
from strucevo.stats.common import stage_rng

AA_BIN_EDGES = (0.0, 0.1, 0.2, 0.35, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 20.0001)
S_BIN_LO = (0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
BLOCK = 20
CORE = 70.0
MIN_N_BIN = 20       # kappa / rho bins with fewer eligible pairs are treated as empty (then filled)
MIN_N_CELL = 5       # lddt SE-model cells with fewer pairs are treated as empty (then filled)


def load_bootsub(root, aln):
    base = Path(root) / "data" / "interim" / "pairs_bootsub"
    dfs, cols = [], []
    for f in sorted(glob.glob(str(base / aln / "shard_[0-9][0-9][0-9][0-9].parquet"))):
        k = Path(f).stem[6:10]
        df = pl.read_parquet(f)
        if (df["split"] != "tune").any():
            raise RuntimeError(f"{f}: boot_sub must contain split == 'tune' pairs only")
        z = np.load(base / aln / f"counts_shard_{k}.npz")
        c = np.load(base / "columns" / aln / f"cols_shard_{k}.npz")
        assert list(z["pair_id"]) == df["pair_id"].to_list() == list(c["pair_id"]), f"order mismatch shard {k}"
        df = df.with_columns(pl.lit(len(cols)).alias("_shard"), pl.Series("_row", np.arange(df.height)))
        dfs.append(df)
        cols.append((dict(z), dict(c)))
    if not dfs:
        raise FileNotFoundError(f"no boot_sub shards for {aln} under {base}")
    return pl.concat(dfs), cols


def pair_columns(cols, shard, row):
    z, c = cols[shard]
    s, e = c["offsets"][row], c["offsets"][row + 1]
    g = {k: c[k][s:e] for k in ("aa_q", "aa_t", "tdi_q", "tdi_t", "plddt_q", "plddt_t", "lddt")}
    core = (g["plddt_q"] >= CORE) & (g["plddt_t"] >= CORE)
    ca = np.where(core & (g["aa_q"] >= 0) & (g["aa_t"] >= 0), g["aa_q"].astype(np.int64) * 20 + g["aa_t"], -1)
    ct = np.where(core & (g["tdi_q"] >= 0) & (g["tdi_t"] >= 0), g["tdi_q"].astype(np.int64) * 20 + g["tdi_t"], -1)
    return ca, ct, g["lddt"].astype(np.float64), int(c["len_a"][row]), z["aa"][row], z["tdi"][row]


def bootstrap_chunk(items, B, rng, dev):
    """items: list of (code_aa, code_tdi, lddt, len_a). Returns counts_aa, counts_tdi (P*B,20,20), S* (P,B)."""
    P = len(items)
    n = np.array([len(it[0]) for it in items])
    Lmax = int(n.max())
    nb = (Lmax + BLOCK - 1) // BLOCK
    ca = np.full((P, Lmax), -1, np.int64); ct = np.full((P, Lmax), -1, np.int64); ll = np.zeros((P, Lmax))
    for i, (a, t, l, _) in enumerate(items):
        ca[i, :n[i]] = a; ct[i, :n[i]] = t; ll[i, :n[i]] = l
    starts = np.floor(rng.random((P, B, nb)) * n[:, None, None]).astype(np.int64)
    T = lambda x: torch.as_tensor(x, device=dev)
    ca, ct, ll, starts, nt = T(ca), T(ct), T(ll), T(starts), T(n)
    j = torch.arange(Lmax, device=dev)
    idx = (starts[:, :, j // BLOCK] + (j % BLOCK)[None, None, :]) % nt[:, None, None]      # (P,B,Lmax)
    valid = j[None, None, :] < nt[:, None, None]
    pb = torch.arange(P * B, device=dev).reshape(P, B, 1)
    out = []
    for codes in (ca, ct):
        g = torch.gather(codes[:, None, :].expand(P, B, Lmax), 2, idx)
        m = valid & (g >= 0)
        flat = (pb * 400 + g)[m]
        out.append(torch.bincount(flat, minlength=P * B * 400).reshape(P * B, 20, 20).cpu().numpy())
    lg = torch.gather(ll[:, None, :].expand(P, B, Lmax), 2, idx)
    S = (lg * valid).sum(2).cpu().numpy() / np.array([it[3] for it in items])[:, None]
    return out[0], out[1], S


def corr_rows(x, y):
    """Per-row Pearson correlation over finite entries of both."""
    m = np.isfinite(x) & np.isfinite(y)
    k = m.sum(1)
    xm = np.where(m, x, 0); ym = np.where(m, y, 0)
    mx = xm.sum(1) / np.maximum(k, 1); my = ym.sum(1) / np.maximum(k, 1)
    dx = np.where(m, x - mx[:, None], 0); dy = np.where(m, y - my[:, None], 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = (dx * dy).sum(1) / np.sqrt((dx ** 2).sum(1) * (dy ** 2).sum(1))
    return np.where(k >= 10, r, np.nan)


def binned(values, b, nbins, agg):
    out, ns = [], []
    for i in range(nbins):
        v = values[(b == i) & np.isfinite(values)]
        ns.append(int(len(v)))
        out.append(float(agg(v)) if len(v) >= MIN_N_BIN else None)
    return out, ns


def zmean(r):
    r = np.clip(r, -0.999999, 0.999999)
    return np.tanh(np.mean(np.arctanh(r)))


def calibrate(root, aln, alpha, B, chunk_cols, dev, log):
    t0 = time.time()
    df, cols = load_bootsub(root, aln)
    df = df.filter(pl.col("status") == "ok")
    s0, s0_src = read_s0(root, aln)
    items, Caa, Ctdi = [], [], []
    n_cnt_mismatch = 0
    S_dump = []
    for sh, row in zip(df["_shard"].to_list(), df["_row"].to_list()):
        ca, ct, ll, la, c0a, c0t = pair_columns(cols, sh, row)
        # consistency: dumped columns reproduce the stored core counts
        if not (np.array_equal(np.bincount(ca[ca >= 0], minlength=400).reshape(20, 20), c0a)
                and np.array_equal(np.bincount(ct[ct >= 0], minlength=400).reshape(20, 20), c0t)):
            n_cnt_mismatch += 1
        items.append((ca, ct, ll, la)); Caa.append(c0a); Ctdi.append(c0t)
        S_dump.append(ll.sum() / la)
    S_dump = np.array(S_dump)
    S = df["lddt_global"].to_numpy()
    lddt_dev = float(np.nanmax(np.abs(S_dump - S))) if len(S) else None
    N = len(items)
    t_load = time.time() - t0
    # point estimates at tuned alpha
    pa = ml_from_counts(np.stack(Caa), "aa", alpha["aa"], dev)
    pt = ml_from_counts(np.stack(Ctdi), "tdi", alpha["tdi"], dev)
    ra = ml_distance(np.stack(Caa), models()["aa"], alpha=alpha["aa"], d_max=D_MAX, device=dev)
    rt = ml_distance(np.stack(Ctdi), models()["tdi"], alpha=alpha["tdi"], d_max=D_MAX, device=dev)
    # bootstrap in length-sorted chunks
    rng = stage_rng("se_calibration", ALNS.index(aln) + 1)
    order = np.argsort([len(it[0]) for it in items], kind="stable")
    dA = np.full((N, B), np.nan); dT = np.full((N, B), np.nan); Sb = np.full((N, B), np.nan)
    t_boot = t_ml = 0.0
    i = 0
    while i < N:
        j = i
        while j < N and (j - i + 1) * B * max(len(items[order[j]][0]), 1) <= chunk_cols:
            j += 1
        j = max(j, i + 1)
        sel = order[i:j]
        t1 = time.time()
        cA, cT, Ss = bootstrap_chunk([items[s] for s in sel], B, rng, dev)
        t2 = time.time()
        for C, key, dst in ((cA, "aa", dA), (cT, "tdi", dT)):
            r = ml_distance(C, models()[key], alpha=alpha[key], d_max=D_MAX, device=dev)
            d = np.where(r.sat_cap, np.nan, r.d)
            dst[sel] = d.reshape(len(sel), B)
        Sb[sel] = Ss
        t_boot += t2 - t1; t_ml += time.time() - t2
        i = j
    capA = np.isnan(dA).mean(1); capT = np.isnan(dT).mean(1)
    sdA = np.nanstd(dA, 1, ddof=1); sdT = np.nanstd(dT, 1, ddof=1)
    la = df["len_a"].to_numpy()
    Db = lddt_D(Sb.ravel(), np.repeat(la, B), s0)["lddt_D"].reshape(N, B)
    sdS = np.std(Sb, 1, ddof=1)
    eligA = np.isfinite(ra.se) & ~ra.saturated & (ra.d > 0) & (capA <= 0.1)
    eligT = np.isfinite(rt.se) & ~rt.saturated & (rt.d > 0) & (capT <= 0.1)
    kA = np.where(eligA, sdA / pa["se"], np.nan)
    kT = np.where(eligT, sdT / pt["se"], np.nan)
    rAT = np.where(eligA & eligT, corr_rows(dA, dT), np.nan)
    rAD = np.where(eligA, corr_rows(dA, Db), np.nan)
    b = bin_index(pa["ml"], AA_BIN_EDGES[:-1])
    nb = len(AA_BIN_EDGES) - 1
    kap_aa, n_kaa = binned(kA, b, nb, np.median)
    kap_tdi, n_ktdi = binned(kT, b, nb, np.median)
    rho_at, n_rat = binned(rAT, b, nb, zmean)
    rho_ad, n_rad = binned(rAD, b, nb, zmean)
    # lddt SE model: sd_S by len bin x S bin
    lb = bin_index(la, LEN_BIN_LO); sb = bin_index(S, S_BIN_LO)
    grid = np.full((len(LEN_BIN_LO), len(S_BIN_LO)), np.nan); ngrid = np.zeros_like(grid, dtype=int)
    for p in range(len(LEN_BIN_LO)):
        for q in range(len(S_BIN_LO)):
            v = sdS[(lb == p) & (sb == q) & np.isfinite(sdS)]
            ngrid[p, q] = len(v)
            if len(v) >= MIN_N_CELL:
                grid[p, q] = np.median(v)
    filled = ~np.isfinite(grid)
    for p in range(len(LEN_BIN_LO)):          # within length bin: nearest S bin
        grid[p] = np.array(fill_nearest(grid[p].tolist())[0], dtype=float)
    for q in range(len(S_BIN_LO)):            # whole empty length bins: nearest length bin
        grid[:, q] = np.array(fill_nearest(grid[:, q].tolist())[0], dtype=float)
    # delta-method check: SD(D*) vs sd_S(model)/(S - S0)
    s0b = np.asarray(s0["median"])[bin_index(la, s0["bin_lo"])]
    se_model = grid[lb, sb] / (S - s0b)
    sdD = np.nanstd(np.where(np.isfinite(Db), Db, np.nan), 1, ddof=1)
    okD = np.isfinite(sdD) & (S - s0b > 0.05) & np.isfinite(se_model)
    fill = lambda v: fill_nearest(v)
    ka_f, ka_flag = fill(kap_aa); kt_f, kt_flag = fill(kap_tdi)
    ra_f, ra_flag = fill(rho_at); rd_f, rd_flag = fill(rho_ad)
    secs = time.time() - t0
    log(f"{aln}: N={N} B={B} load {t_load:.1f}s boot {t_boot:.1f}s ml {t_ml:.1f}s total {secs:.1f}s; "
        f"count mismatches {n_cnt_mismatch}, max|S_dump-S| {lddt_dev}")
    log(f"  kappa_aa {np.round(np.array(kap_aa, dtype=float), 2)} n {n_kaa}")
    log(f"  kappa_tdi {np.round(np.array(kap_tdi, dtype=float), 2)} n {n_ktdi}")
    log(f"  rho_aa_tdi {np.round(np.array(rho_at, dtype=float), 2)}  rho_aa_lddt {np.round(np.array(rho_ad, dtype=float), 2)}")
    yaml_part = {
        "kappa": {"bin_edges": list(AA_BIN_EDGES), "aa_ml": ka_f, "tdi_ml": kt_f},
        "rho": {"bin_edges": list(AA_BIN_EDGES), "aa_tdi": ra_f, "aa_lddt": rd_f},
        "lddt_se_model": {"len_bin_lo": list(LEN_BIN_LO), "s_bin_lo": list(S_BIN_LO), "sd_S": grid.tolist()}}
    det = {"n_pairs_ok": N, "B": B, "block": BLOCK, "alpha": alpha, "s0_source": s0_src,
           "kappa_raw": {"aa_ml": kap_aa, "tdi_ml": kap_tdi, "n_aa": n_kaa, "n_tdi": n_ktdi,
                         "filled_aa": ka_flag, "filled_tdi": kt_flag,
                         "pooled_median_aa": float(np.nanmedian(kA)) if np.isfinite(kA).any() else None,
                         "pooled_median_tdi": float(np.nanmedian(kT)) if np.isfinite(kT).any() else None},
           "rho_raw": {"aa_tdi": rho_at, "aa_lddt": rho_ad, "n_aa_tdi": n_rat, "n_aa_lddt": n_rad,
                       "filled_aa_tdi": ra_flag, "filled_aa_lddt": rd_flag},
           "lddt_se_model_n": ngrid.tolist(), "lddt_se_model_filled": filled.tolist(),
           "lddt_se_model_min_n_cell": MIN_N_CELL, "kappa_rho_min_n_bin": MIN_N_BIN,
           "delta_method_check": {"n": int(okD.sum()),
                                  "median_ratio_sdD_boot_over_model": float(np.median(sdD[okD] / se_model[okD]))
                                  if okD.any() else None},
           "frac_pairs_capped_reps_gt10pct": {"aa": float((capA > 0.1).mean()), "tdi": float((capT > 0.1).mean())},
           "checks": {"count_mismatches": n_cnt_mismatch, "max_abs_S_dump_minus_table": lddt_dev},
           "seconds": {"load": round(t_load, 1), "bootstrap_counts": round(t_boot, 1), "ml": round(t_ml, 1),
                       "total": round(secs, 1)}}
    return yaml_part, det


def iid_sim_check(alpha_aa, dev, n_pairs=600, B=100, seed=3):
    """Machinery check: columns simulated i.i.d. under LG+G4(alpha_aa) -> kappa should be ~1."""
    from scipy.linalg import expm
    from strucevo.distances import gamma_rates
    m = models()["aa"]
    Q = m.Q(); r, _ = gamma_rates(alpha_aa, 4)
    rng = np.random.default_rng(seed)
    items = []
    for _ in range(n_pairs):
        d = rng.uniform(0.1, 2.0); L = int(rng.integers(200, 401))
        cat = rng.integers(0, 4, L); a = rng.choice(20, L, p=m.pi)
        cum = [np.cumsum(expm(Q * rc * d), 1) for rc in r]
        b = np.array([min(np.searchsorted(cum[c][x], rng.random() * cum[c][x][-1]), 19) for c, x in zip(cat, a)])
        code = a * 20 + b
        items.append((code, code, np.ones(L), L))
    cA, _, _ = bootstrap_chunk(items, B, rng, dev)
    C0 = np.stack([np.bincount(it[0], minlength=400).reshape(20, 20) for it in items])
    r0 = ml_distance(C0, m, alpha=alpha_aa, d_max=D_MAX, device=dev)
    rb = ml_distance(cA, m, alpha=alpha_aa, d_max=D_MAX, device=dev).d.reshape(n_pairs, B)
    k = np.std(rb, 1, ddof=1) / r0.se
    return {"n_pairs": n_pairs, "B": B, "alpha": alpha_aa, "kappa_median": float(np.nanmedian(k)),
            "kappa_mean": float(np.nanmean(k))}


def main():
    ap = add_root_arg(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter))
    ap.add_argument("--alns", nargs="+", default=list(ALNS))
    ap.add_argument("--B", type=int, default=100)
    ap.add_argument("--alpha-aa", type=float)
    ap.add_argument("--alpha-tdi", type=float)
    ap.add_argument("--chunk-cols", type=int, default=20_000_000, help="P*B*Lmax per GPU chunk")
    ap.add_argument("--device", default=None)
    a = ap.parse_args()
    dev = a.device or ("cuda" if torch.cuda.is_available() else "cpu")
    alpha, src = read_alpha(a.root, a.alpha_aa, a.alpha_tdi)
    res = {"step": "se_calibration", "split": "tune", "split_md5": split_md5(a.root), "alpha_used": alpha,
           "alpha_source": src, "kappa": {}, "rho": {}, "lddt_se_model": {}, "details": {}}
    for aln in a.alns:
        y, det = calibrate(a.root, aln, alpha, a.B, a.chunk_cols, dev, lambda s: print(s, flush=True))
        for k, v in y.items():
            res[k][aln] = v
        res["details"][aln] = det
    res["details"]["iid_sim_check"] = iid_sim_check(min(alpha["aa"], 50.0), dev)
    print("iid simulation check (kappa should be ~1):", res["details"]["iid_sim_check"], flush=True)
    res["provenance"] = provenance()
    print(write_json(res, out_json(a.root, "se_calibration")))


if __name__ == "__main__":
    main()
