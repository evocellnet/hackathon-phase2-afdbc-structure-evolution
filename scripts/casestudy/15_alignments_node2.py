"""X6 node 2: extra alignment sources for the alignment-sensitivity of branch lengths.

famsa : FAMSA 2 (default) AA MSA of the 400 sampled sequences; 3Di projected residue-by-residue (as 06_alignments.py)
fmg   : FoldMason easy-msa with --guide-tree = the (outgroup-rooted, binary) AA ML tree. The default FoldMason actin
        MSA (06) isolates 67-69 FtsA/MreB rows in their own columns (all-gap after the 50% filter); guided by the
        ML tree it does not (0 rows < 10% occupancy). Note: the guide tree is the tree being fixed (mild circularity,
        acceptable because the topology is fixed anyway).
Column filter occupancy >= 50% (same as 06). Appends rows to results/casestudy/alignment_stats.tsv.
"""
import numpy as np, pandas as pd
R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
CS = f"{R}/data/interim/casestudy"


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


stats = pd.read_csv(f"{R}/results/casestudy/alignment_stats.tsv", sep="\t")
stats = stats[~stats.source.isin(["famsa", "fmg"])]
new = []
for sf in ("actin", "tubulin"):
    raw = rf(f"{CS}/{sf}_sample.fasta")
    tdi = rf(f"{CS}/{sf}_3di.fasta")
    fa = {k: v.upper() for k, v in rf(f"{CS}/famsa/{sf}_famsa.fasta").items()}
    td = {}
    for k, s in fa.items():
        assert s.replace("-", "") == raw[k]
        it = iter(tdi[k]); td[k] = "".join("-" if c == "-" else next(it) for c in s)
    fmA = rf(f"{CS}/foldmason/guide/{sf}_aa.fa"); fmT = rf(f"{CS}/foldmason/guide/{sf}_3di.fa")
    for k in fmA:
        assert fmA[k].replace("-", "") == raw[k], k
    for src, (A, T) in {"famsa": (fa, td), "fmg": (fmA, fmT)}.items():
        keys = list(raw)
        M = np.array([list(A[k]) for k in keys])
        keep = (M != "-").mean(0) >= 0.5
        At = {k: "".join(np.array(list(A[k]))[keep]) for k in keys}
        Tt = {k: "".join(np.array(list(T[k]))[keep]) for k in keys}
        assert all((a == "-") == (b == "-") for k in keys for a, b in zip(At[k], Tt[k]))
        wf(At, f"{CS}/aln/{sf}_{src}_aa.fasta"); wf(Tt, f"{CS}/aln/{sf}_{src}_3di.fasta")
        wf(A, f"{CS}/aln/{sf}_{src}_aa.full.fasta"); wf(T, f"{CS}/aln/{sf}_{src}_3di.full.fasta")
        new.append(dict(superfamily=sf, source=src, n_seq=len(keys), cols_full=M.shape[1], cols_kept=int(keep.sum()),
                        mean_gap_frac_kept=round(np.mean([s.count("-") / len(s) for s in At.values()]), 3),
                        median_seq_len=int(np.median([len(raw[k]) for k in keys]))))
out = pd.concat([stats, pd.DataFrame(new)])
out.to_csv(f"{R}/results/casestudy/alignment_stats.tsv", sep="\t", index=False)
print(out.to_string())
