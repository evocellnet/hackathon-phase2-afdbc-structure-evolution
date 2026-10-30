"""X7 step 6: alignments of the sampled kinase domains and column filtering.

Sources (all on the same <= 500 cropped domains):
  linsi     : MAFFT L-INS-i (--localpair --maxiterate 1000, 8 threads) AA MSA  [primary; fixed-topology reference]
  famsa     : FAMSA 2 (default) AA MSA                                         [sensitivity]
  foldmason : FoldMason easy-msa on the cropped domain models (AA + 3Di on shared columns) [sensitivity]
              (optionally foldmason_guided with --guide-tree = AA ML tree, if the default isolates rows)
For linsi/famsa the 3Di letters are projected residue-by-residue onto the AA columns.
Column filter: occupancy >= 50 %, identical for AA and 3Di. Rows with < 10 % occupancy after filtering are reported.
Usage: python 06_alignments.py run|filter [src ...]
"""
import os, subprocess, sys
import numpy as np, pandas as pd

R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
K = f"{R}/data/interim/kinase"
os.makedirs(f"{K}/aln", exist_ok=True)


def rf(p):
    d, k = {}, None
    for line in open(p):
        line = line.strip()
        if line.startswith(">"):
            k = line[1:].split()[0].replace(".pdb", ""); d[k] = []
        elif k:
            d[k].append(line)
    return {k: "".join(v) for k, v in d.items()}


def wf(d, p):
    with open(p, "w") as fo:
        for k, v in d.items():
            fo.write(f">{k}\n{v}\n")


mode = sys.argv[1]
srcs = sys.argv[2:] or ["linsi", "famsa", "foldmason"]
if mode == "run":
    fa = f"{K}/pkl_sample.fasta"
    if "linsi" in srcs and not os.path.exists(f"{K}/aln/linsi_raw.fasta"):
        subprocess.run(f"mafft --localpair --maxiterate 1000 --thread 8 --quiet {fa} > {K}/aln/linsi_raw.fasta",
                       shell=True, check=True)
    if "famsa" in srcs and not os.path.exists(f"{K}/aln/famsa_raw.fasta"):
        subprocess.run(f"famsa -t 8 {fa} {K}/aln/famsa_raw.fasta", shell=True, check=True)
    if "foldmason" in srcs and not os.path.exists(f"{K}/foldmason/pkl_aa.fa"):
        os.makedirs(f"{K}/foldmason", exist_ok=True)
        subprocess.run(f"foldmason easy-msa {K}/pdb {K}/foldmason/pkl {K}/foldmason/tmp --threads 8 -v 1",
                       shell=True, check=True)
    if "foldmason_guided" in srcs and not os.path.exists(f"{K}/foldmason/guided_aa.fa"):
        subprocess.run(f"foldmason easy-msa {K}/pdb {K}/foldmason/guided {K}/foldmason/tmpg --threads 8 -v 1 "
                       f"--guide-tree {K}/iqtree/pkl_guide.nwk", shell=True, check=True)
    sys.exit()

raw = rf(f"{K}/pkl_sample.fasta")
tdi = rf(f"{K}/pkl_3di.fasta")
pdbaa = rf(f"{K}/pkl_aa_from_pdb.fasta")
assert set(raw) == set(tdi) == set(pdbaa)
bad = [k for k in raw if raw[k] != pdbaa[k]]
assert not bad, f"{len(bad)} AA/PDB mismatches e.g. {bad[:3]}"
statf = f"{R}/results/kinase/alignment_stats.tsv"
stats = pd.read_csv(statf, sep="\t") if os.path.exists(statf) else pd.DataFrame()
if len(stats):
    stats = stats[~stats.source.isin(srcs)]
new = []
for src in srcs:
    if src in ("linsi", "famsa"):
        A = {k: v.upper() for k, v in rf(f"{K}/aln/{src}_raw.fasta").items()}
        T = {}
        for k, s in A.items():
            assert s.replace("-", "") == raw[k], k
            it = iter(tdi[k]); T[k] = "".join("-" if c == "-" else next(it) for c in s)
    else:
        pre = "pkl" if src == "foldmason" else "guided"
        A = rf(f"{K}/foldmason/{pre}_aa.fa"); T = rf(f"{K}/foldmason/{pre}_3di.fa")
        for k in A:
            assert A[k].replace("-", "") == raw[k], k
    keys = list(raw)
    M = np.array([list(A[k]) for k in keys])
    keep = (M != "-").mean(0) >= 0.5
    At = {k: "".join(np.array(list(A[k]))[keep]) for k in keys}
    Tt = {k: "".join(np.array(list(T[k]))[keep]) for k in keys}
    assert all((a == "-") == (b == "-") for k in keys for a, b in zip(At[k], Tt[k]))
    rowocc = {k: 1 - At[k].count("-") / len(At[k]) for k in keys}
    lowrows = [k for k, v in rowocc.items() if v < 0.10]
    wf(At, f"{K}/aln/{src}_aa.fasta"); wf(Tt, f"{K}/aln/{src}_3di.fasta")
    wf(A, f"{K}/aln/{src}_aa.full.fasta"); wf(T, f"{K}/aln/{src}_3di.full.fasta")
    new.append(dict(source=src, n_seq=len(keys), cols_full=M.shape[1], cols_kept=int(keep.sum()),
                    mean_gap_frac_kept=round(np.mean([s.count("-") / len(s) for s in At.values()]), 3),
                    n_rows_occ_lt10pct=len(lowrows), median_dom_len=int(np.median([len(raw[k]) for k in keys]))))
    print(src, new[-1], lowrows[:10])
out = pd.concat([stats, pd.DataFrame(new)])
out.to_csv(statf, sep="\t", index=False)
print(out.to_string())
