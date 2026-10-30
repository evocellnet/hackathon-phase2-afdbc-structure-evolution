"""X6 step 3b: column filtering and 3Di projection.

For each superfamily and each alignment source:
  linsi     : MAFFT L-INS-i AA MSA (fixed-topology reference); 3Di letters are projected onto the AA columns
              residue-by-residue (AFDB v4 model sequence == UniProt sequence, verified here).
  foldmason : FoldMason structure-guided MSA (AA and 3Di outputs share columns by construction).
Column filter: keep columns with occupancy >= 50 % (documented, applied identically to AA and 3Di).
Outputs data/interim/casestudy/aln/{sf}_{src}_{aa,3di}.fasta and results/casestudy/alignment_stats.tsv
"""
import os
import numpy as np, pandas as pd

R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
CS = f"{R}/data/interim/casestudy"
OCC = 0.5


def read_fa(p):
    d, k = {}, None
    for line in open(p):
        line = line.strip()
        if line.startswith(">"):
            k = line[1:].split()[0]
            d[k] = []
        elif k:
            d[k].append(line)
    return {k: "".join(v) for k, v in d.items()}


def write_fa(d, p):
    with open(p, "w") as fo:
        for k, v in d.items():
            fo.write(f">{k}\n{v}\n")


os.makedirs(f"{CS}/aln", exist_ok=True)
stats = []
for sf in ["actin", "tubulin"]:
    raw = read_fa(f"{CS}/{sf}_sample.fasta")
    tdi = {k.split()[0]: v for k, v in read_fa(f"{CS}/{sf}_3di.fasta").items()}
    pdbaa = read_fa(f"{CS}/{sf}_aa_from_pdb.fasta")
    assert set(raw) == set(tdi) == set(pdbaa)
    nbad = sum(raw[k] != pdbaa[k] for k in raw)
    assert nbad == 0, f"{nbad} AA/PDB mismatches"
    assert all(len(raw[k]) == len(tdi[k]) for k in raw)
    # --- MAFFT: project 3Di
    aa = read_fa(f"{CS}/{sf}_linsi.fasta")
    aa = {k: v.upper() for k, v in aa.items()}
    td = {}
    for k, s in aa.items():
        assert s.replace("-", "") == raw[k]
        it = iter(tdi[k])
        td[k] = "".join("-" if c == "-" else next(it) for c in s)
    srcs = {"linsi": (aa, td)}
    # --- FoldMason
    fm_aa = read_fa(f"{CS}/foldmason/{sf}_aa.fa")
    fm_td = read_fa(f"{CS}/foldmason/{sf}_3di.fa")
    fm_aa = {k.replace(".pdb", ""): v for k, v in fm_aa.items()}
    fm_td = {k.replace(".pdb", ""): v for k, v in fm_td.items()}
    for k in fm_aa:
        assert fm_aa[k].replace("-", "") == raw[k], k
    srcs["foldmason"] = (fm_aa, fm_td)
    for src, (A, T) in srcs.items():
        keys = list(raw)
        M = np.array([list(A[k]) for k in keys])
        occ = (M != "-").mean(0)
        keep = occ >= OCC
        At = {k: "".join(np.array(list(A[k]))[keep]) for k in keys}
        Tt = {k: "".join(np.array(list(T[k]))[keep]) for k in keys}
        write_fa(At, f"{CS}/aln/{sf}_{src}_aa.fasta")
        write_fa(Tt, f"{CS}/aln/{sf}_{src}_3di.fasta")
        write_fa(A, f"{CS}/aln/{sf}_{src}_aa.full.fasta")
        write_fa(T, f"{CS}/aln/{sf}_{src}_3di.full.fasta")
        gapfrac = np.mean([s.count("-") / len(s) for s in At.values()])
        stats.append(dict(superfamily=sf, source=src, n_seq=len(keys), cols_full=M.shape[1],
                          cols_kept=int(keep.sum()), mean_gap_frac_kept=round(gapfrac, 3),
                          median_seq_len=int(np.median([len(raw[k]) for k in keys]))))
pd.DataFrame(stats).to_csv(f"{R}/results/casestudy/alignment_stats.tsv", sep="\t", index=False)
print(pd.DataFrame(stats).to_string())
