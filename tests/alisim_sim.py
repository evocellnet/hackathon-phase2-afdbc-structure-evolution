"""Regenerable AliSim simulations + IQ-TREE pairwise ML distances for validating strucevo.ml_distance.

For each (model, length, d): AliSim simulates `n_rep` two-taxon alignments on the tree (A:d,B:0) under
MODEL+G4{alpha}. All replicate pairs are then stacked into one 2*n_rep-sequence alignment and IQ-TREE is run
with the same fixed model (-t BIONJ -n 0, chunks of 20 pairs) to obtain its pairwise ML distances (.mldist); only within-pair entries are used.

Outputs (gitignored, regenerable): data/interim/tests/alisim/<model>_L<len>_d<d>/
    sim_<i>.fa        AliSim replicates
    iq_<k>.log        IQ-TREE logs (chunks of CHUNK stacked pairs, -t BIONJ -n 0 -> .mldist)
    summary.npz       true d, IQ-TREE distances, sequences
Usage: python tests/alisim_sim.py [--jobs 24]
"""
from __future__ import annotations

import argparse
import glob
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "data", "interim", "tests", "alisim")
GLOBAL_SEED = 20261030
DS = [0.05, 0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0]
LENGTHS = [300, 1000]
MODELS = ["LG", "WAG"]
ALPHA = 1.0
N_REP = 200
CHUNK = 20  # pairs per IQ-TREE run
IQTREE = os.environ.get("IQTREE", "iqtree3")


def cfg_dir(model, L, d):
    return os.path.join(OUT, f"{model}_L{L}_d{d:g}")


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def read_mldist(path):
    with open(path) as fh:
        n = int(fh.readline().split()[0])
        names, rows = [], []
        for _ in range(n):
            parts = fh.readline().split()
            names.append(parts[0])
            rows.append([float(x) for x in parts[1:]])
    return names, np.array(rows)


def run_one(args):
    model, L, d, seed = args
    wd = cfg_dir(model, L, d)
    if os.path.exists(os.path.join(wd, "summary.npz")):
        return wd
    os.makedirs(wd, exist_ok=True)
    with open(os.path.join(wd, "tree.nwk"), "w") as fh:
        fh.write(f"(A:{d},B:0);\n")
    mstr = f"{model}+G4{{{ALPHA}}}"
    subprocess.run([IQTREE, "--alisim", os.path.join(wd, "sim"), "-t", os.path.join(wd, "tree.nwk"), "-m", mstr,
                    "--length", str(L), "--num-alignments", str(N_REP), "--seed", str(seed), "-af", "fasta",
                    "-redo", "-nt", "1"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    A, B = [], []
    for i in range(1, N_REP + 1):
        s = read_fasta(os.path.join(wd, f"sim_{i}.fa"))
        A.append(s["A"]); B.append(s["B"])
    # IQ-TREE needs >= 3 sequences: stack CHUNK pairs per run; -t BIONJ makes it compute + write ML distances
    # under the fixed model (no free parameters with LG/WAG+G4{alpha}); -n 0 skips tree search.
    iq = np.full(N_REP, np.nan)
    for c0 in range(0, N_REP, CHUNK):
        ids = range(c0, min(c0 + CHUNK, N_REP))
        fa = os.path.join(wd, f"pairs_{c0}.fa")
        with open(fa, "w") as fh:
            for i in ids:
                fh.write(f">A{i}\n{A[i]}\n>B{i}\n{B[i]}\n")
        pre = os.path.join(wd, f"iq_{c0}")
        subprocess.run([IQTREE, "-s", fa, "-m", mstr, "-t", "BIONJ", "-n", "0", "-nt", "1", "-keep-ident",
                        "-seed", str(seed), "--prefix", pre, "-redo", "-quiet"],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        names, M = read_mldist(pre + ".mldist")
        pos = {n: k for k, n in enumerate(names)}
        for i in ids:
            iq[i] = M[pos[f"A{i}"], pos[f"B{i}"]]
        for f in glob.glob(pre + ".*") + [fa]:
            if not f.endswith(".log"):
                os.remove(f)
    np.savez(os.path.join(wd, "summary.npz"), d_true=d, length=L, model=model, iq=iq,
             A=np.array(A), B=np.array(B), seed=seed)
    for f in glob.glob(os.path.join(wd, "sim_*.fa")):
        os.remove(f)
    return wd


def all_configs():
    out = []
    for mi, model in enumerate(MODELS):
        for li, L in enumerate(LENGTHS):
            for di, d in enumerate(DS):
                out.append((model, L, d, GLOBAL_SEED + 1000 * mi + 100 * li + di))
    return out


def ensure(jobs=16):
    cfgs = all_configs()
    todo = [c for c in cfgs if not os.path.exists(os.path.join(cfg_dir(*c[:3]), "summary.npz"))]
    if todo:
        import multiprocessing as mp
        with ProcessPoolExecutor(jobs, mp_context=mp.get_context("spawn")) as ex:
            list(ex.map(run_one, todo))
    return [cfg_dir(m, L, d) for m, L, d, _ in cfgs]


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=16)
    a = ap.parse_args()
    for wd in ensure(a.jobs):
        print(wd, file=sys.stderr)
