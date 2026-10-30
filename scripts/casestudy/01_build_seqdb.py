"""Build the sequence DB for the actin/tubulin case study (X6).
All panel proteins (180 euk + 100 prok) whose AFDB v4 model sequence equals the
UniProt 2026_03 sequence (md5 match). pLDDT filter is applied later so that
the membership table can report how many were lost to it.
Header: >acc|leaf
"""
import pandas as pd, glob, os
R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
out = f"{R}/data/interim/casestudy/panel_all_v4match.fasta"
keep = set()
tabs = []
for f in ["protein_table.parquet", "protein_table_prok.parquet"]:
    d = pd.read_parquet(f"{R}/data/processed/panel/{f}",
                        columns=["acc", "leaf", "superkingdom", "seq_match", "mean_plddt", "len_v4"])
    tabs.append(d)
t = pd.concat(tabs)
t = t[t.seq_match]
keep = set(t.acc + "|" + t.leaf)
n = 0
with open(out, "w") as fo:
    for fa in sorted(glob.glob(f"{R}/data/interim/structures/v4_fasta/*.fasta")):
        leaf = os.path.basename(fa)[:-6]
        hdr = None
        for line in open(fa):
            if line.startswith(">"):
                acc = line[1:].split()[0]
                k = f"{acc}|{leaf}"
                hdr = k if k in keep else None
                if hdr:
                    fo.write(f">{hdr}\n"); n += 1
            elif hdr:
                fo.write(line)
print("written", n, "of", len(keep))
