#!/usr/bin/env python3
"""Sanity checks for the panel data (writes results/panel/sanity_checks_<name>.json).

1. 20 random proteins with a fetched AFDB v4 model (seed 20261030): live UniProt REST sequence vs local 2026
   proteome FASTA sequence vs FoldComp model sequence.
2. Whole-table agreement between the two independent v4-sequence checks (TED md5 vs md5 of the decoded
   FoldComp sequence).
3. Fetch completeness: proteins with a v4 lookup entry vs proteins decoded from the local FoldComp DB.
4. Per-species counts (min / median / max) and species below the 80% exact-match gate.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys
import time
import urllib.request

import polars as pl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import ROOT  # noqa: E402

SEED = 20261030


def uniprot_seq(acc, tries=6):
    req = urllib.request.Request(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", headers={"User-Agent": "Mozilla/5.0"})
    for i in range(tries):  # UniProt REST returns sporadic HTTP 500s
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                txt = r.read().decode()
            return "".join(txt.split("\n")[1:]).strip(), txt.split("\n")[0]
        except Exception:  # noqa: BLE001
            if i == tries - 1:
                raise
            time.sleep(3 * (i + 1))


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "panel_v4"
    table = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "data/processed/panel/protein_table.parquet")
    tab = pl.read_parquet(table)
    st = pl.read_parquet(os.path.join(ROOT, "data/interim/structures", name + "_proteins.parquet"))
    out = {"name": name, "table": os.path.relpath(table, ROOT)}
    # 1. random 20
    smp = tab.filter(pl.col("fetched")).sample(20, seed=SEED).join(st.select("acc", "seq_v4"), on="acc")
    rows = []
    for r in smp.iter_rows(named=True):
        try:
            s_live, hdr = uniprot_seq(r["acc"])
        except Exception as e:  # noqa: BLE001
            s_live, hdr = None, str(e)
        md5_live = hashlib.md5(s_live.encode()).hexdigest() if s_live else None
        rows.append({"acc": r["acc"], "leaf": r["leaf"], "len_2026": r["len_2026"], "len_v4": r["len_v4"],
                     "live_eq_local2026": md5_live == r["md5_2026"], "foldcomp_eq_local2026": r["seq_match"],
                     "foldcomp_eq_live": s_live == r["seq_v4"], "uniprot_header": hdr[:100]})
    out["random20"] = rows
    out["random20_foldcomp_eq_live"] = sum(x["foldcomp_eq_live"] for x in rows)
    out["random20_live_eq_local2026"] = sum(x["live_eq_local2026"] for x in rows)
    # 2. TED vs FoldComp
    f = tab.filter(pl.col("fetched"))
    out["n_fetched"] = f.height
    out["n_seq_match_foldcomp"] = int(f["seq_match"].sum())
    out["n_seq_match_ted"] = int(f["seq_match_ted"].fill_null(False).sum())
    out["n_ted_foldcomp_disagree"] = int((f["seq_match"] != f["seq_match_ted"].fill_null(False)).sum())
    out["mismatch_rate_among_fetched"] = round(1 - out["n_seq_match_foldcomp"] / max(1, f.height), 5)
    mm = f.filter(~pl.col("seq_match"))
    out["mismatch_len_diff_zero_frac"] = round(float((mm["len_2026"] == mm["len_v4"]).mean()), 4) if mm.height else None
    # 3. completeness
    out["n_canonical"] = tab.height
    out["n_has_v4"] = int(tab["has_v4"].sum())
    out["n_has_v4_not_fetched"] = int((tab["has_v4"] & ~tab["fetched"]).sum())
    out["n_no_v4"] = int((~tab["has_v4"]).sum())
    out["n_no_v4_len_gt2700"] = int(((~tab["has_v4"]) & (tab["len_2026"] > 2700)).sum())
    out["n_len_gt2700_total"] = int((tab["len_2026"] > 2700).sum())
    out["frac_proteins_seq_match"] = round(out["n_seq_match_foldcomp"] / tab.height, 4)
    # 4. per species
    sp = tab.group_by("leaf").agg(pl.len().alias("n"), pl.col("seq_match").mean().alias("frac_match"),
                                  pl.col("mean_plddt").mean().alias("plddt"))
    out["n_species"] = sp.height
    out["per_species_n"] = {k: float(v) for k, v in zip(["min", "median", "max"], [sp["n"].min(), sp["n"].median(), sp["n"].max()])}
    out["per_species_frac_match"] = {k: round(float(v), 4) for k, v in zip(["min", "median", "max"],
                                     [sp["frac_match"].min(), sp["frac_match"].median(), sp["frac_match"].max()])}
    out["species_below_0.80"] = sp.filter(pl.col("frac_match") < 0.80).sort("frac_match").rows()
    out["mean_plddt_overall"] = round(float(tab["mean_plddt"].mean()), 2)
    os.makedirs(os.path.join(ROOT, "results/panel"), exist_ok=True)
    with open(os.path.join(ROOT, "results/panel", f"sanity_checks_{name}.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "random20"}, indent=1))


if __name__ == "__main__":
    main()
