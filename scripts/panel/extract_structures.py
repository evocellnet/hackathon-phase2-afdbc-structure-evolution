#!/usr/bin/env python3
"""Decode every entry of the local panel FoldComp DB and build the per-protein table.

Outputs
  data/interim/structures/<name>_proteins.parquet : acc, entry, fc_key, seq_v4, len_v4, md5_v4, mean_plddt,
                                                    frac_plddt_ge70, frac_plddt_ge90
  data/interim/structures/<name>_plddt.npz       : per-residue pLDDT (float16; FoldComp stores B-factors
                                                    lossily, ~0.1 precision): acc[], offsets[], plddt[]
  data/interim/structures/v4_fasta/<leaf>.fasta  : v4 model sequences per species (header = accession)
  data/processed/panel/protein_table.parquet     : one row per canonical protein of the panel(s):
        acc, upid, leaf, species name, superkingdom, len_2026, md5_2026, has_v4, len_v4, md5_v4,
        seq_match (md5_2026 == md5_v4 from the FoldComp model), seq_match_ted, mean_plddt, frac_plddt_ge70
  data/interim/panel/structure_coverage.tsv       : per species counts (canonical, has_v4, fetched, seq_match)
"""
from __future__ import annotations

import argparse
import hashlib
import mmap
import multiprocessing as mp
import os
import sys

import numpy as np
import polars as pl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import ROOT  # noqa: E402

_MM = None


def _init(path):
    global _MM
    fh = open(path, "rb")
    _MM = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)


def _decode(chunk):
    import foldcomp
    out = []
    for key, off, ln, entry in chunk:
        b = _MM[off:off + ln]
        try:
            d = foldcomp.get_data(b)
        except Exception as e:  # noqa: BLE001
            out.append((key, entry, None, None, str(e)))
            continue
        seq = d["residues"]
        pl_ = np.asarray(d["b_factors"], dtype=np.float32)
        out.append((key, entry, seq, pl_, None))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--name", default="panel_v4")
    ap.add_argument("--panel", nargs="+", required=True, help="panel TSV(s) defining the protein_table rows")
    ap.add_argument("--table", default=os.path.join(ROOT, "data/processed/panel/protein_table.parquet"))
    ap.add_argument("--procs", type=int, default=46)
    a = ap.parse_args()
    sdir = os.path.join(ROOT, "data/interim/structures")
    db = os.path.join(sdir, a.name)
    idx = pl.read_csv(db + ".index", separator="\t", has_header=False, new_columns=["fc_key", "off", "len"])
    lk = pl.read_csv(db + ".lookup", separator="\t", has_header=False, new_columns=["fc_key", "entry", "x"])
    idx = idx.join(lk.select("fc_key", "entry"), on="fc_key")
    rows = list(idx.iter_rows())
    chunks = [rows[i:i + 2000] for i in range(0, len(rows), 2000)]
    recs, errs = [], []
    accs, offs, vals = [], [0], []
    with mp.Pool(a.procs, initializer=_init, initargs=(db,)) as pool:
        for res in pool.imap(_decode, chunks):
            for key, entry, seq, pll, err in res:
                acc = entry.split("-")[1]
                if err:
                    errs.append((acc, err))
                    continue
                recs.append((acc, entry, key, seq, len(seq), hashlib.md5(seq.encode()).hexdigest(),
                             float(pll.mean()), float((pll >= 70).mean()), float((pll >= 90).mean())))
                accs.append(acc)
                vals.append(pll.astype(np.float16))
                offs.append(offs[-1] + len(pll))
    print(f"decoded {len(recs)} entries; {len(errs)} decode errors {errs[:5]}", flush=True)
    st = pl.DataFrame(recs, schema=["acc", "entry", "fc_key", "seq_v4", "len_v4", "md5_v4", "mean_plddt",
                                    "frac_plddt_ge70", "frac_plddt_ge90"], orient="row")
    st.write_parquet(os.path.join(sdir, a.name + "_proteins.parquet"))
    np.savez_compressed(os.path.join(sdir, a.name + "_plddt.npz"), acc=np.array(accs), offsets=np.array(offs, dtype=np.int64),
                        plddt=np.concatenate(vals) if vals else np.zeros(0, np.float16))
    # per-protein table
    leaves = pl.concat([pl.read_csv(p, separator="\t").select("leaf", "name", "superkingdom") for p in a.panel])
    prot = pl.read_parquet(os.path.join(ROOT, "data/interim/panel/proteins_2026.parquet")).join(leaves, on="leaf")
    tab = (prot.join(st.select("acc", "len_v4", "md5_v4", "mean_plddt", "frac_plddt_ge70", "frac_plddt_ge90")
                     .rename({"len_v4": "len_v4_model"}), on="acc", how="left")
           .with_columns(pl.col("md5_v4").is_not_null().alias("fetched"),
                         (pl.col("md5_v4") == pl.col("md5_2026")).fill_null(False).alias("seq_match")))
    tab = tab.select("acc", "upid", "leaf", "name", "superkingdom", "db", "len_2026", "md5_2026", "has_v4", "fetched",
                     pl.col("len_v4_model").alias("len_v4"), "md5_v4", "seq_match", "seq_match_ted", "mean_plddt",
                     "frac_plddt_ge70", "frac_plddt_ge90")
    os.makedirs(os.path.dirname(a.table), exist_ok=True)
    tab.write_parquet(a.table)
    cov = (tab.group_by("leaf", "name").agg(
        pl.len().alias("n_canonical"), pl.col("has_v4").sum().alias("n_has_v4"), pl.col("fetched").sum().alias("n_fetched"),
        pl.col("seq_match").sum().alias("n_seq_match"), (pl.col("fetched") & ~pl.col("seq_match")).sum().alias("n_mismatch"),
        (~pl.col("has_v4")).sum().alias("n_no_v4"), ((~pl.col("has_v4")) & (pl.col("len_2026") > 2700)).sum().alias("n_no_v4_gt2700"),
        (pl.col("seq_match") != pl.col("seq_match_ted").fill_null(False)).sum().alias("n_ted_vs_foldcomp_disagree"),
        pl.col("mean_plddt").mean().round(2).alias("mean_plddt"))
        .with_columns((pl.col("n_seq_match") / pl.col("n_canonical")).round(4).alias("frac_seq_match")).sort("leaf"))
    cov.write_csv(os.path.join(ROOT, "data/interim/panel", f"structure_coverage_{a.name}.tsv"), separator="\t")
    # v4 model-sequence FASTAs per species
    fdir = os.path.join(sdir, "v4_fasta")
    os.makedirs(fdir, exist_ok=True)
    j = tab.filter(pl.col("fetched")).select("leaf", "acc").join(st.select("acc", "seq_v4"), on="acc")
    for (leaf,), g in j.group_by(["leaf"]):
        with open(os.path.join(fdir, f"{leaf}.fasta"), "w") as fo:
            for acc, s in g.select("acc", "seq_v4").sort("acc").iter_rows():
                fo.write(f">{acc}\n{s}\n")
    print(cov.select(pl.col("n_canonical", "n_has_v4", "n_fetched", "n_seq_match", "n_mismatch", "n_no_v4", "n_no_v4_gt2700",
                            "n_ted_vs_foldcomp_disagree").sum()))


if __name__ == "__main__":
    main()
