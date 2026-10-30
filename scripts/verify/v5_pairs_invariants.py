#!/usr/bin/env python3
"""V5 invariants of the per-pair tables (ANALYSIS_PLAN §7; the age check lives with the panel/timetree code).

  * rows out == rows in, per alignment (vs the input pair list); pair_id unique; gene_a is the shorter protein
  * self pairs (gene_a == gene_b): TM (tm_short; Foldseek alignments) = 1 and lddt_global = 1 (|Δ| <= 1e-3);
    aa_ml = tdi_ml = 0 where n_core > 0; aa_ml_seq = 0
  * identical sequences (md5_a == md5_b, includes self pairs): aa_ml_seq = 0 (pass criterion, all tables).
    Core aa_ml of the structural alignment is reported (n > 0) but not a pass criterion: TM-align can
    misregister identical sequences whose two AF models differ in conformation (found in the euk20 demo).
  * md5 identity: md5(UniProt FASTA sequence) == md5(AFDB v4 FoldComp sequence) per store protein (rate reported;
    the D1 filter acts on this at HOG construction)
  * human-mouse 1:1 orthologs: median AA identity (aa_pid_seq) within [0.80, 0.90]

  python scripts/verify/v5_pairs_invariants.py --pairs-in P.tsv --out data/interim/pairs_demo/pairs \
      --store S --fasta-glob 'data/interim/ancestra/euk20/*.fasta.gz' --json results/verify/V5_pairs.json
"""
from __future__ import annotations

import argparse
import glob
import gzip
import hashlib
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "src"))

import numpy as np  # noqa: E402
import polars as pl  # noqa: E402

from strucevo.pairs.store import Store  # noqa: E402

HUMAN, MOUSE = 9606, 10090


def fasta_md5(paths, wanted):
    out = {}
    for p in paths:
        with gzip.open(p, "rt") as fh:
            acc, seq = None, []
            for line in fh:
                if line.startswith(">"):
                    if acc in wanted:
                        out[acc] = hashlib.md5("".join(seq).encode()).hexdigest()
                    acc = line[1:].split("|")[1] if "|" in line else line[1:].split()[0]
                    seq = []
                else:
                    seq.append(line.strip())
            if acc in wanted:
                out[acc] = hashlib.md5("".join(seq).encode()).hexdigest()
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pairs-in", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--store", required=True)
    ap.add_argument("--fasta-glob")
    ap.add_argument("--json", required=True)
    a = ap.parse_args()
    store = Store(a.store)
    md5 = dict(zip(store.meta["acc"].to_list(), store.meta["md5"].to_list()))
    n_in = pl.read_csv(a.pairs_in, separator="\t", infer_schema_length=0).height
    res = {"n_pairs_in": n_in}
    ok = True
    for aln in ("seqonly", "tmalign", "tdiaa"):
        fs = glob.glob(os.path.join(a.out, aln, "shard_[0-9][0-9][0-9][0-9].parquet"))
        if not fs:
            continue
        d = pl.read_parquet(fs)
        r = {"rows_out": d.height, "rows_equal": d.height == n_in, "pair_id_unique": d["pair_id"].n_unique() == d.height,
             "gene_a_shorter": bool((d["len_a"] <= d["len_b"]).all()),
             "status": dict(d.group_by("status").len().iter_rows())}
        s = d.filter(pl.col("gene_a") == pl.col("gene_b"))
        r["n_self"] = s.height
        r["self_lddt_global_min"] = float(s["lddt_global"].min())
        r["self_tm_short_min"] = float(s["tm_short"].min()) if aln != "seqonly" else None
        sc = s.filter(pl.col("n_core") > 0)
        r["self_aa_ml_max"] = float(sc["aa_ml"].max()); r["self_tdi_ml_max"] = float(sc["tdi_ml"].max())
        r["self_aa_ml_seq_max"] = float(s["aa_ml_seq"].max())
        r["self_ok"] = bool(r["self_lddt_global_min"] >= 0.999 and r["self_aa_ml_max"] == 0 and r["self_tdi_ml_max"] == 0
                            and r["self_aa_ml_seq_max"] == 0
                            and (aln == "seqonly" or abs(r["self_tm_short_min"] - 1) <= 1e-3))
        ma = np.array([md5.get(g) for g in d["gene_a"].to_list()], dtype=object)
        mb = np.array([md5.get(g) for g in d["gene_b"].to_list()], dtype=object)
        ident = d.filter(pl.Series(ma == mb))
        r["n_identical_seq"] = ident.height
        r["n_identical_nonself"] = int((ident["gene_a"] != ident["gene_b"]).sum())
        r["identical_aa_ml_seq_max"] = float(ident["aa_ml_seq"].max())
        r["identical_aa_ml_max"] = float(ident.filter(pl.col("n_core") > 0)["aa_ml"].max())
        r["identical_n_core_aa_ml_gt0"] = int((ident.filter(pl.col("n_core") > 0)["aa_ml"] > 0).sum())
        r["identical_ok"] = bool(r["identical_aa_ml_seq_max"] == 0)
        r["identical_warn"] = r["identical_n_core_aa_ml_gt0"] > 0
        hm = d.filter((pl.col("one_to_one") == 1) & (((pl.col("species_a") == HUMAN) & (pl.col("species_b") == MOUSE)) |
                                                     ((pl.col("species_a") == MOUSE) & (pl.col("species_b") == HUMAN))))
        r["human_mouse_1to1"] = {"n": hm.height, "median_aa_pid_seq": float(hm["aa_pid_seq"].median()),
                                 "median_aa_ml_seq": float(hm["aa_ml_seq"].median()),
                                 "median_aa_ml_core_ok": float(hm.filter(pl.col("status") == "ok")["aa_ml"].median()),
                                 "median_aa_p_core_ok": float(hm.filter(pl.col("status") == "ok")["aa_p"].median())}
        r["human_mouse_ok"] = bool(0.80 <= r["human_mouse_1to1"]["median_aa_pid_seq"] <= 0.90)
        ok &= r["rows_equal"] and r["pair_id_unique"] and r["gene_a_shorter"] and r["self_ok"] and r["identical_ok"]
        res[aln] = r
    if a.fasta_glob:
        accs = set(store.meta["acc"].to_list())
        fm = fasta_md5(sorted(glob.glob(a.fasta_glob)), accs)
        n_found = len(fm)
        n_eq = sum(fm[k] == md5[k] for k in fm)
        res["md5_identity"] = {"n_store": len(accs), "n_in_fasta": n_found, "n_equal": n_eq,
                               "mismatch_rate": (n_found - n_eq) / max(n_found, 1)}
    res["V5_pass"] = bool(ok and all(res[x]["human_mouse_ok"] for x in ("seqonly",) if x in res))
    os.makedirs(os.path.dirname(os.path.abspath(a.json)), exist_ok=True)
    json.dump(res, open(a.json, "w"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
