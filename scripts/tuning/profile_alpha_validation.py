#!/usr/bin/env python3
"""Validate strucevo.distances.profile_alpha (pooled Gamma-shape estimate from two-sequence count matrices).

Two simulation sources at known alpha, both with discrete Gamma (4 category means, as IQ-TREE/AliSim +G4):
  (1) AliSim (iqtree3 --alisim; independent implementation): tree (A:d,B:0), MODEL+G4{alpha}, --length L,
      n_rep alignments per (model, alpha, d, L); models LG and Q.3Di.AF (PAML file, ARND order -> AA_ORDER).
  (2) in-house strucevo.distances.simulate_counts with d ~ U(0.1, 2), L ~ U{200..400}.
For each (source, model, alpha) the pooled profile likelihood (per-pair d free; d_max = 20) is maximised with
the same fit_alpha() used by tune_alpha.py; recovered alpha, profile 95% CI and the relative bias are reported,
also for d sub-ranges (0.1-0.5, 0.5-2). A check passes if true alpha is inside the CI or |bias| <= 5%.

  python scripts/tuning/profile_alpha_validation.py [--root R] [--n-rep 60] [--jobs 2]
Output: ROOT/results/tuning/profile_alpha_validation.json, ROOT/figures/tuning/profile_alpha_validation.png
Simulations (regenerable) in ROOT/data/interim/tuning/alisim/.  Fragment has no tuning.yaml keys (validation only).
"""
from __future__ import annotations

import argparse
import os
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from se_common import PROJECT, add_root_arg, fit_alpha, models, np, out_json, provenance, write_json

from strucevo.distances import AA_ORDER, ml_distance, pair_counts, simulate_counts

DS = (0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0)
LS = (200, 300, 400)
ALPHAS = (0.5, 1.0, 2.0)
SEED = 20261030


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def alisim_one(args):
    wd, mstr, d, L, n_rep, seed = args
    done = os.path.join(wd, "pairs.npz")
    if os.path.exists(done):
        return done
    os.makedirs(wd, exist_ok=True)
    tree = os.path.join(wd, "tree.nwk")
    open(tree, "w").write(f"(A:{d},B:0);\n")
    subprocess.run(["iqtree3", "--alisim", os.path.join(wd, "sim"), "-t", tree, "-m", mstr, "--length", str(L),
                    "--num-alignments", str(n_rep), "--seed", str(seed), "-af", "fasta", "-redo", "-nt", "1",
                    "--seqtype", "AA"],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    A, B = [], []
    for i in range(1, n_rep + 1):
        f = os.path.join(wd, f"sim_{i}.fa")
        s = read_fasta(f)
        A.append(s["A"]); B.append(s["B"])
        os.remove(f)
    C = pair_counts(A, B, AA_ORDER)
    np.savez_compressed(done, counts=C, d=np.full(n_rep, d), L=np.full(n_rep, L))
    return done


def summarise(C, d, model, alpha_true):
    out = {}
    for name, sel in (("all", np.ones(len(d), bool)), ("d0.1-0.5", d <= 0.5), ("d0.5-2", d > 0.5)):
        r = fit_alpha(C[sel], model)
        lo, hi = r["ci95"]
        bias = r["alpha"] / alpha_true - 1
        out[name] = dict(alpha_hat=r["alpha"], ci95=[lo, hi], rel_bias=bias, n_pairs=int(sel.sum()),
                         covered=bool(lo <= alpha_true <= hi), ok=bool(lo <= alpha_true <= hi or abs(bias) <= 0.05))
    # distance recovery at the true alpha (sanity of the engine)
    res = ml_distance(C, model, alpha=alpha_true, d_max=20)
    out["d_rel_bias_median_at_true_alpha"] = float(np.median(res.d / d - 1))
    res_hat = ml_distance(C, model, alpha=out["all"]["alpha_hat"], d_max=20)
    out["d_rel_bias_median_at_alpha_hat"] = float(np.median(res_hat.d / d - 1))
    return out


def main():
    ap = add_root_arg(argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter))
    ap.add_argument("--n-rep", type=int, default=60)
    ap.add_argument("--n-inhouse", type=int, default=5000)
    ap.add_argument("--jobs", type=int, default=2)
    a = ap.parse_args()
    t0 = time.time()
    root = Path(a.root)
    simdir = root / "data" / "interim" / "tuning" / "alisim"
    M = models()
    q3di = str(PROJECT / "data" / "raw" / "models" / "Q.3Di.AF")
    mspec = {"LG": ("LG", M["aa"]), "Q.3Di.AF": (q3di, M["tdi"])}
    res = {"step": "profile_alpha_validation", "design": {"d": DS, "L": LS, "alpha": ALPHAS, "n_rep": a.n_rep,
                                                           "n_inhouse": a.n_inhouse, "ncat": 4, "d_max": 20},
           "results": {}}
    jobs = []
    k = 0
    for mname, (mpath, _) in mspec.items():
        for al in ALPHAS:
            for d in DS:
                for L in LS:
                    k += 1
                    wd = str(simdir / f"{mname}_a{al:g}_d{d:g}_L{L}")
                    jobs.append((wd, f"{mpath}+G4{{{al}}}", d, L, a.n_rep, SEED + k))
    t1 = time.time()
    with ThreadPoolExecutor(a.jobs) as ex:
        list(ex.map(alisim_one, jobs))
    res["alisim_seconds"] = round(time.time() - t1, 1)
    rng = np.random.default_rng(SEED)
    for mname, (_, model) in mspec.items():
        for al in ALPHAS:
            Cs, ds = [], []
            for d in DS:
                for L in LS:
                    z = np.load(simdir / f"{mname}_a{al:g}_d{d:g}_L{L}" / "pairs.npz")
                    Cs.append(z["counts"]); ds.append(z["d"])
            C, d = np.concatenate(Cs), np.concatenate(ds)
            key = f"{mname}/alpha{al:g}"
            r_ali = summarise(C, d, model, al)
            dd = rng.uniform(0.1, 2.0, a.n_inhouse)
            LL = rng.integers(200, 401, a.n_inhouse)
            Ci = simulate_counts(model, dd, LL, alpha=al, rng=rng)
            r_in = summarise(Ci, dd, model, al)
            res["results"][key] = {"alisim": r_ali, "inhouse": r_in}
            print(key, "alisim", round(r_ali["all"]["alpha_hat"], 3), [round(x, 3) for x in r_ali["all"]["ci95"]],
                  "| inhouse", round(r_in["all"]["alpha_hat"], 3), [round(x, 3) for x in r_in["all"]["ci95"]],
                  flush=True)
    allok = [v[s][sub]["ok"] for v in res["results"].values() for s in ("alisim", "inhouse") for sub in ("all",)]
    res["verdict"] = {"pooled_all_ok": bool(all(allok)),
                      "max_abs_rel_bias_pooled": float(max(abs(v[s]["all"]["rel_bias"]) for v in res["results"].values()
                                                         for s in ("alisim", "inhouse")))}
    res["provenance"] = provenance()
    res["seconds"] = round(time.time() - t0, 1)
    print(write_json(res, out_json(a.root, "profile_alpha_validation")))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 2, figsize=(9, 3.6), sharey=True)
        for ax, mname in zip(axs, mspec):
            for j, src in enumerate(("alisim", "inhouse")):
                for sub, mk in (("all", "o"), ("d0.1-0.5", "v"), ("d0.5-2", "^")):
                    x = np.array(ALPHAS) * (1 + 0.04 * (j - 0.5) + 0.015 * ("all d0.1-0.5 d0.5-2".split().index(sub) - 1))
                    y = [res["results"][f"{mname}/alpha{al:g}"][src][sub]["alpha_hat"] for al in ALPHAS]
                    lo = [res["results"][f"{mname}/alpha{al:g}"][src][sub]["ci95"][0] for al in ALPHAS]
                    hi = [res["results"][f"{mname}/alpha{al:g}"][src][sub]["ci95"][1] for al in ALPHAS]
                    ax.errorbar(x, y, yerr=[np.subtract(y, lo), np.subtract(hi, y)], fmt=mk, ms=4,
                                color=("C0", "C1")[j], label=f"{src} {sub}", alpha=0.85)
            ax.plot([0.4, 2.5], [0.4, 2.5], "k:", lw=1)
            ax.set_xscale("log"); ax.set_yscale("log")
            ax.set_title(mname); ax.set_xlabel("true alpha")
        axs[0].set_ylabel("recovered pooled alpha (95% profile CI)")
        axs[1].legend(fontsize=6)
        fig.tight_layout()
        fp = root / "figures" / "tuning" / "profile_alpha_validation.png"
        fp.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(fp, dpi=150)
        print(fp)
    except Exception as e:  # figure optional
        print("figure failed:", e)


if __name__ == "__main__":
    main()
