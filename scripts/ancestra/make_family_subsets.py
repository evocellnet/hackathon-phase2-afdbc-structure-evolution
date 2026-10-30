#!/usr/bin/env python3
"""Make species-subsampled FASTA files of one protein family from a UniProt TSV export
(columns Entry, Organism (ID), ..., Sequence). Header: ">ACC OX=<taxid>".
Usage: python scripts/ancestra/make_family_subsets.py FAMILY.tsv.gz OUTPREFIX N1 [N2 ...]   (N = species; 0 = all)
Species are sampled uniformly at random (seed 11); nested subsets (smaller ones are subsets of larger)."""
import random, sys
import pandas as pd
tsv, pre, *ns = sys.argv[1:]
df = pd.read_csv(tsv, sep="\t")
sp = sorted(df["Organism (ID)"].unique())
random.Random(11).shuffle(sp)
for n in map(int, ns):
    keep = set(sp if n == 0 else sp[:n])
    sub = df[df["Organism (ID)"].isin(keep)]
    out = f"{pre}_sp{n if n else len(sp)}.fasta"
    with open(out, "w") as fh:
        for acc, t, s in sub[["Entry", "Organism (ID)", "Sequence"]].itertuples(index=False):
            fh.write(f">{acc} OX={t}\n{s}\n")
    print(out, len(keep), len(sub))
