#!/usr/bin/env python3
"""Per-protein gene-model flags for the panel (ANALYSIS_PLAN §1.6 step 1).

Output: data/interim/hogs/gene_flags.parquet with columns
  acc, frag_uniprot (bool; '(Fragment)' in the UniProt 2026_03 FASTA description),
  cluflag (int; AFDB Clusters v3 file 5 cluFlag: 1 AFDB50 member, 2 AFDB-cluster member,
           3 fragment, 4 singleton; 0 = accession absent from file 5),
  afdb_repid (str)
Run once (~5 min; single core + zcat).
"""
import gzip, subprocess, sys, os, pandas as pd, numpy as np
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
PT = ROOT / "data/processed/panel/protein_table.parquet"
F5 = ROOT / "data/raw/afdb_clusters/v3/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz"
OUT = ROOT / "data/interim/hogs/gene_flags.parquet"

def main(table=PT, out=OUT, proteome_dir=ROOT / "data/raw/proteomes"):
    pt = pd.read_parquet(table, columns=["acc", "leaf"])
    accs = set(pt.acc)
    # 1) UniProt fragment flag from FASTA descriptions
    frag = set()
    for leaf in sorted(pt.leaf.unique()):
        fa = proteome_dir / f"{leaf}.fasta.gz"
        if not fa.exists():
            sys.exit(f"missing proteome {fa}")
        with gzip.open(fa, "rt") as fh:
            for line in fh:
                if line[0] == ">":
                    acc = line[1:].split("|")[1] if "|" in line else line[1:].split()[0]
                    if "(Fragment)" in line:
                        frag.add(acc)
    print(f"UniProt '(Fragment)' entries in panel: {len(frag & accs)}", flush=True)
    # 2) AFDB cluFlag (stream file 5 through awk keyed on the panel accessions)
    tmp = out.with_suffix(".accs.txt")
    tmp.write_text("\n".join(sorted(accs)) + "\n")
    cmd = (f"zcat {F5} | awk -F'\\t' 'FNR==NR{{a[$1]=1;next}} ($2 in a){{print $2\"\\t\"$3\"\\t\"$1}}' {tmp} -")
    res = subprocess.run(["bash", "-c", cmd], capture_output=True, text=True, check=True)
    rows = [l.split("\t") for l in res.stdout.splitlines()]
    cf = pd.DataFrame(rows, columns=["acc", "cluflag", "afdb_repid"]).drop_duplicates("acc")
    cf["cluflag"] = cf.cluflag.astype(int)
    tmp.unlink()
    df = pd.DataFrame({"acc": sorted(accs)})
    df["frag_uniprot"] = df.acc.isin(frag)
    df = df.merge(cf, on="acc", how="left")
    df["cluflag"] = df.cluflag.fillna(0).astype(int)
    df.to_parquet(out, index=False)
    print(df.frag_uniprot.sum(), df.cluflag.value_counts().to_dict(), flush=True)

if __name__ == "__main__":
    main()
