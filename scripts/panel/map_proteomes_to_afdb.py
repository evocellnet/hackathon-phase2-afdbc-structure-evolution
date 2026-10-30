#!/usr/bin/env python3
"""Map downloaded canonical proteomes (UniProt current release) to AFDB v4 (FoldComp afdb_uniprot_v4).

For every protein: accession, proteome, 2026 sequence length + md5, whether an AFDB v4 model exists
(FoldComp lookup: entry AF-<acc>-F1-model_v4), and whether the v4 model's sequence equals the current
sequence (md5 from the TED per-chain table, which lists md5(sequence) for all 214M AFDB v4 models; verified
on a FoldComp-decoded entry).  Writes
  data/interim/panel/proteins_2026.parquet       (one row per canonical protein of every downloaded proteome)
  data/interim/panel/proteome_v4_coverage.tsv    (per proteome counts; read by select_species.py)
Optionally (--ancestra-panel TSV --ancestra-out DIR) writes Ancestra input FASTAs (2026 canonical sequences,
header = accession, file = <UPID>_<taxid>.fasta) for the species of a panel.
"""
from __future__ import annotations

import argparse
import glob
import gzip
import hashlib
import os
import re
import sys

import polars as pl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import ROOT  # noqa: E402

PROT_DIR = os.path.join(ROOT, "data/raw/proteomes")
INTERIM = os.path.join(ROOT, "data/interim/panel")


def read_fasta(path):
    acc, seq, hdr = None, [], None
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if acc is not None:
                    yield acc, hdr, "".join(seq)
                hdr = line[1:].rstrip()
                tok = hdr.split()[0]
                acc = tok.split("|")[1] if tok.count("|") >= 2 else tok
                seq = []
            else:
                seq.append(line.strip())
    if acc is not None:
        yield acc, hdr, "".join(seq)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ancestra-panel", nargs="*", default=[])
    ap.add_argument("--ancestra-out")
    ap.add_argument("--fasta-only", action="store_true", help="only write the Ancestra FASTAs")
    a = ap.parse_args()
    if not a.fasta_only:
        map_all()
    if a.ancestra_out:
        write_ancestra(a.ancestra_panel, a.ancestra_out)


def map_all():
    rows = []
    files = sorted(glob.glob(os.path.join(PROT_DIR, "UP*_*.fasta.gz")))
    for f in files:
        leaf = os.path.basename(f).replace(".fasta.gz", "")
        up, tx = leaf.split("_")
        for acc, hdr, s in read_fasta(f):
            m = re.search(r" OX=(\d+)", hdr)
            rows.append((acc, up, leaf, int(tx), int(m.group(1)) if m else None, hdr.split("|")[0],
                         len(s), hashlib.md5(s.upper().encode()).hexdigest()))
    df = pl.DataFrame(rows, schema=["acc", "upid", "leaf", "proteome_taxid", "ox", "db", "len_2026", "md5_2026"],
                      orient="row")
    print(f"{len(files)} proteomes, {df.height} proteins", flush=True)
    accs = df.select("acc").unique()
    lk = (pl.scan_parquet(os.path.join(INTERIM, "afdb_v4_lookup.parquet"))
          .join(accs.lazy(), on="acc", how="semi").select("acc", "fc_key", "entry").collect())
    ted = (pl.scan_parquet(os.path.join(INTERIM, "afdb_v4_ted_md5.parquet"))
           .join(accs.lazy(), on="acc", how="semi").select("acc", pl.col("md5").alias("md5_v4_ted"),
                                                           pl.col("nres").alias("len_v4")).collect())
    df = (df.join(lk, on="acc", how="left").join(ted.unique("acc"), on="acc", how="left")
            .with_columns(pl.col("fc_key").is_not_null().alias("has_v4"))
            .with_columns((pl.col("has_v4") & (pl.col("md5_v4_ted") == pl.col("md5_2026"))).alias("seq_match_ted")))
    df.write_parquet(os.path.join(INTERIM, "proteins_2026.parquet"))
    cov = (df.group_by("upid", "leaf").agg(
        pl.len().alias("n_canonical"), pl.col("has_v4").sum().alias("n_v4"),
        pl.col("seq_match_ted").sum().alias("n_v4_exact"),
        (pl.col("has_v4") & pl.col("md5_v4_ted").is_null()).sum().alias("n_v4_no_ted_md5"),
        (~pl.col("has_v4")).sum().alias("n_missing"),
        ((~pl.col("has_v4")) & (pl.col("len_2026") > 2700)).sum().alias("n_missing_gt2700"),
        (pl.col("has_v4") & ~pl.col("seq_match_ted").fill_null(False)).sum().alias("n_v4_mismatch"),
        (pl.col("db") == "sp").sum().alias("n_swissprot"))
        .with_columns((pl.col("n_v4_exact") / pl.col("n_canonical")).round(4).alias("frac_v4_exact"),
                      (pl.col("n_v4") / pl.col("n_canonical")).round(4).alias("frac_v4"))
        .sort("upid"))
    cov.write_csv(os.path.join(INTERIM, "proteome_v4_coverage.tsv"), separator="\t")
    print(cov["frac_v4_exact"].describe())


def write_ancestra(panels, outdir):
    if True:
        os.makedirs(outdir, exist_ok=True)
        leaves = set()
        for p in panels:
            leaves |= set(pl.read_csv(p, separator="\t")["leaf"].to_list())
        for leaf in sorted(leaves):
            src = os.path.join(PROT_DIR, leaf + ".fasta.gz")
            dst = os.path.join(outdir, leaf + ".fasta")
            with open(dst + ".tmp", "w") as fo:
                for acc, hdr, s in read_fasta(src):
                    fo.write(f">{acc}\n{s}\n")
            os.replace(dst + ".tmp", dst)
        print(f"wrote {len(leaves)} Ancestra FASTAs to {outdir}")


if __name__ == "__main__":
    main()
