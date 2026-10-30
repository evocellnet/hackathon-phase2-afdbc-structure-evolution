"""X6 node 2, step 1a: IQ-TREE branch lengths on a fixed topology, alignment sensitivity and column bootstrap.

For each (superfamily, alignment source, tree):
  full  : iqtree3 -te <tree> -m LG+G4 on the AA columns and -m Q.3Di.AF+G4 on the 3Di letters of the SAME columns
  boot  : B replicates resampling alignment columns with replacement (the SAME column indices for AA and 3Di, so
          the AA/3Di covariance from shared columns is kept); model and Gamma alpha fixed at the full-data estimate
          (only branch lengths re-optimised).
Trees: 'ml'   = the binary AA ML tree (data/interim/casestudy/iqtree/{sf}_linsi_aa_ml.treefile; UFBoot<95 branches
                are nuisance parameters, not interpreted)
       'cons' = the constrained ML tree (13_constraint_trees.py; data/interim/casestudy/constraint/{sf}_cons.treefile)
Outputs data/interim/casestudy/iqfix/{sf}_{src}_{tree}_{aa,3di}_full.treefile and boot/{...}_b{k}.treefile
Usage: python 16_iqtree_fixed.py --sf actin --src linsi famsa fmg --tree ml --B 100 --workers 8
"""
import argparse, os, re, subprocess, shutil
from concurrent.futures import ThreadPoolExecutor
import numpy as np

R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
CS = f"{R}/data/interim/casestudy"
Q3DI = f"{R}/data/raw/models/Q.3Di.AF"
OUT = f"{CS}/iqfix"
TREES = {"ml": f"{CS}/iqtree/{{sf}}_linsi_aa_ml.treefile", "cons": f"{CS}/constraint/{{sf}}_cons.treefile"}
ap = argparse.ArgumentParser()
ap.add_argument("--sf", nargs="+", default=["actin", "tubulin"])
ap.add_argument("--src", nargs="+", default=["linsi"])
ap.add_argument("--tree", default="ml")
ap.add_argument("--B", type=int, default=100)
ap.add_argument("--workers", type=int, default=8)
ap.add_argument("--model_aa", default="LG+G4")
ap.add_argument("--tag", default="")
ap.add_argument("--full_threads", type=int, default=4)
a = ap.parse_args()
os.makedirs(f"{OUT}/boot", exist_ok=True)


def rf(p):
    d, k = {}, None
    for line in open(p):
        line = line.strip()
        if line.startswith(">"):
            k = line[1:].split()[0]; d[k] = []
        elif k:
            d[k].append(line)
    return {k: "".join(v) for k, v in d.items()}


def run(aln, tree, model, prefix, threads=1):
    if os.path.exists(prefix + ".treefile"):
        return prefix + ".treefile"
    cmd = ["iqtree3", "-s", aln, "-te", tree, "-m", model, "-T", str(threads), "--seed", "20261030",
           "--prefix", prefix, "-redo", "--quiet"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(prefix + ".treefile"):
        raise RuntimeError(" ".join(cmd) + "\n" + r.stdout[-2000:] + r.stderr[-2000:])
    for ext in (".ckp.gz", ".log", ".model.gz", ".uniqueseq.phy", ".iqtree"):
        if "/boot/" in prefix and os.path.exists(prefix + ext):
            os.remove(prefix + ext)
    return prefix + ".treefile"


def alpha_of(prefix):
    txt = open(prefix + ".iqtree").read()
    return float(re.search(r"Gamma shape alpha:\s+([0-9.]+)", txt).group(1))


rng = np.random.default_rng(20261030)
jobs = []
for sf in a.sf:
    tree = TREES[a.tree].format(sf=sf)
    for src in a.src:
        tag = f"{sf}_{src}_{a.tree}{a.tag}"
        A = rf(f"{CS}/aln/{sf}_{src}_aa.fasta"); T = rf(f"{CS}/aln/{sf}_{src}_3di.fasta")
        keys = list(A); L = len(A[keys[0]])
        m_aa = a.model_aa
        m_td = f"{Q3DI}+G4"
        pa = run(f"{CS}/aln/{sf}_{src}_aa.fasta", tree, m_aa, f"{OUT}/{tag}_aa_full", threads=a.full_threads)
        pt = run(f"{CS}/aln/{sf}_{src}_3di.fasta", tree, m_td, f"{OUT}/{tag}_3di_full", threads=a.full_threads)
        al_aa, al_td = alpha_of(f"{OUT}/{tag}_aa_full"), alpha_of(f"{OUT}/{tag}_3di_full")
        print(tag, "alpha", al_aa, al_td, flush=True)
        fixed_aa = re.sub(r"G4$", f"G4{{{al_aa}}}", m_aa)
        fixed_td = f"{Q3DI}+G4{{{al_td}}}"
        Aa = np.array([list(A[k]) for k in keys]); Ta = np.array([list(T[k]) for k in keys])
        for b in range(a.B):
            cols = rng.integers(0, L, L)
            for X, nm, mod in ((Aa, "aa", fixed_aa), (Ta, "3di", fixed_td)):
                p = f"{OUT}/boot/{tag}_{nm}_b{b}"
                fa = p + ".fasta"
                if not os.path.exists(p + ".treefile"):
                    with open(fa, "w") as fo:
                        for k, row in zip(keys, X[:, cols]):
                            fo.write(f">{k}\n{''.join(row)}\n")
                jobs.append((fa, tree, mod, p))


def job(j):
    fa, tree, mod, p = j
    out = run(fa, tree, mod, p)
    if os.path.exists(fa):
        os.remove(fa)
    return out


with ThreadPoolExecutor(a.workers) as ex:
    for k, _ in enumerate(ex.map(job, jobs)):
        if k % 50 == 0:
            print(f"boot {k}/{len(jobs)}", flush=True)
print("done")
