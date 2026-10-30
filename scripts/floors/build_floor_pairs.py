#!/usr/bin/env python3
"""Build the floor_random and floor_ident pair lists (ANALYSIS_PLAN §1.2/§1.4), eukaryote panel only.

Candidate proteins: AFDB v4 model fetched, sequence md5 == 2026 sequence (seq_match), mean pLDDT >= 70, len 50..2700.
floor_random: per length bin of the SHORTER protein (<100 (>=50), 100-200, 200-300, 300-500, 500-800, >=800)
  * different species (species_taxid), different AFDB-v3 structural cluster (file 5 repId),
    no shared Pfam unit (unit = clan if the family has one, else the family), both proteins carry >= 1 Pfam family,
  * partner length in [len_a, 2*len_a] (HOG-like length ratio), every protein used <= MAXUSE times.
floor_ident: one random pair of distinct accessions per md5 group (identical AFDB v4 sequence), up to N per length bin.
Seeds: SeedSequence(20261030).spawn(16)[3] (config/seeds.yaml, stage floor_random_pairs).
"""
import argparse, gzip, os, sys, json
import numpy as np, polars as pl

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BINS = [(50, 100), (100, 200), (200, 300), (300, 500), (500, 800), (800, 2701)]
BIN_NAMES = ["<100", "100-200", "200-300", "300-500", "500-800", ">=800"]
TARGET = [10000, 20000, 20000, 20000, 18000, 12000]  # 100k total
MAXUSE = 2


def bin_of(L):
    for i, (lo, hi) in enumerate(BINS):
        if lo <= L < hi:
            return i
    return -1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=os.path.join(ROOT, "data/interim/floors"))
    ap.add_argument("--n-ident-per-bin", type=int, default=3000)
    a = ap.parse_args()
    ss = np.random.SeedSequence(20261030).spawn(16)[3]
    rng_r, rng_i = [np.random.default_rng(s) for s in ss.spawn(2)]
    F = os.path.join(ROOT, "data/interim/floors")

    t = pl.read_parquet(os.path.join(ROOT, "data/processed/panel/protein_table.parquet"))
    sp = pl.read_csv(os.path.join(ROOT, "data/processed/panel/species_panel.tsv"), separator="\t").select("leaf", "taxid", "species_taxid")
    e = (t.filter((pl.col("superkingdom") == "Eukaryota") & pl.col("fetched") & pl.col("seq_match")
                  & (pl.col("mean_plddt") >= 70) & (pl.col("len_v4") <= 2700) & (pl.col("len_v4") >= 50))
         .join(sp, on="leaf", how="left")
         .select("acc", "leaf", "taxid", "species_taxid", pl.col("len_v4").cast(pl.Int64).alias("len"), "md5_v4", "mean_plddt"))
    assert e["taxid"].null_count() == 0
    # clusters
    cl = pl.read_parquet(os.path.join(F, "cand_cluster5.parquet")).rename({"rep": "cluster", "flag": "cluflag"})
    e = e.join(cl, on="acc", how="left")
    assert e["cluster"].null_count() == 0
    # Pfam
    clans = pl.read_csv(os.path.join(ROOT, "data/raw/pfam/Pfam-A.clans.tsv.gz"), separator="\t", has_header=False,
                        new_columns=["pf", "clan", "cn", "fn", "desc"], infer_schema_length=0)
    unit = {r[0]: (r[1] if r[1] else r[0]) for r in clans.select("pf", "clan").iter_rows()}
    h = pl.read_csv(os.path.join(F, "pfam_hits.tsv"), separator="\t", has_header=False,
                    new_columns=["pacc", "md5", "pf", "s", "en"], infer_schema_length=0).select("pacc", "md5", "pf")
    by_acc = h.select(pl.col("pacc").alias("acc"), "pf").unique()
    by_md5 = h.select(pl.col("md5").alias("md5_v4"), "pf").unique()
    pa = e.select("acc", "md5_v4").join(by_acc, on="acc", how="inner").select("acc", "pf")
    pm = e.select("acc", "md5_v4").join(by_md5, on="md5_v4", how="inner").select("acc", "pf")
    pfam = pl.concat([pa, pm]).unique()
    pfam = pfam.with_columns(pl.col("pf").map_elements(lambda x: unit.get(x, x), return_dtype=pl.Utf8).alias("unit"))
    pu = pfam.group_by("acc").agg(pl.col("unit").unique().alias("units"), pl.col("pf").unique().alias("fams"))
    e = e.join(pu, on="acc", how="left")
    e = e.with_columns(pl.col("units").is_not_null().alias("annot"))
    stats = {"n_candidates": e.height, "n_annotated": int(e["annot"].sum())}
    stats["frac_annotated_by_bin"] = {}
    for i, (lo, hi) in enumerate(BINS):
        s = e.filter((pl.col("len") >= lo) & (pl.col("len") < hi))
        stats["frac_annotated_by_bin"][BIN_NAMES[i]] = [s.height, float(s["annot"].mean())]
    print(json.dumps(stats), flush=True)

    # ---------------- floor_random ----------------
    p = e.filter(pl.col("annot")).sort(["len", "acc"])
    L = p["len"].to_numpy()
    sid = p["species_taxid"].to_numpy()
    cid = p["cluster"].to_numpy()
    unit_ids = {}
    U = [frozenset(unit_ids.setdefault(u, len(unit_ids)) for u in us) for us in p["units"].to_list()]
    acc = p["acc"].to_numpy()
    leaf_tax = p["taxid"].to_numpy()
    n = len(L)
    use = np.zeros(n, np.int8)
    seen = set()
    rows = []
    rej = {"species": 0, "cluster": 0, "pfam": 0, "dup": 0, "use": 0}
    for b, (lo, hi) in enumerate(BINS):
        i0, i1 = np.searchsorted(L, lo), np.searchsorted(L, hi)
        got = 0
        tries = 0
        while got < TARGET[b] and tries < 60:
            tries += 1
            ia = rng_r.integers(i0, i1, size=int((TARGET[b] - got) * 1.6) + 100)
            j0 = ia  # partner length >= len_a
            j1 = np.searchsorted(L, 2 * L[ia], side="right")
            ib = j0 + (rng_r.random(len(ia)) * (j1 - j0)).astype(np.int64)
            for x, y in zip(ia, ib):
                if got >= TARGET[b]:
                    break
                if x == y:
                    continue
                if use[x] >= MAXUSE or use[y] >= MAXUSE:
                    rej["use"] += 1; continue
                if sid[x] == sid[y]:
                    rej["species"] += 1; continue
                if cid[x] == cid[y]:
                    rej["cluster"] += 1; continue
                if U[x] & U[y]:
                    rej["pfam"] += 1; continue
                k = (min(x, y), max(x, y))
                if k in seen:
                    rej["dup"] += 1; continue
                seen.add(k)
                use[x] += 1; use[y] += 1
                rows.append((x, y)); got += 1
        print(BIN_NAMES[b], got, "tries", tries, flush=True)
        assert got >= 2000, (b, got)
    rr = np.array(rows)
    x, y = rr[:, 0], rr[:, 1]
    # orient: a shorter (ties lexicographic by acc); p sorted by (len, acc) so index order == orientation
    lo_i, hi_i = np.minimum(x, y), np.maximum(x, y)
    fr = pl.DataFrame({
        "gene_a": acc[lo_i], "gene_b": acc[hi_i],
        "pair_set": "floor_random", "hog": "NA", "split": "NA",
        "species_a": leaf_tax[lo_i], "species_b": leaf_tax[hi_i],
        "relation": "random", "len_a_fc": L[lo_i], "len_b_fc": L[hi_i],
        "bin": [BIN_NAMES[bin_of(v)] for v in L[lo_i]],
    })
    fr.write_csv(os.path.join(F, "floor_random_pairs.tsv"), separator="\t")
    stats["random_rejections"] = rej
    stats["random_pairs_by_bin"] = {k: int(v) for k, v in fr.group_by("bin").len().iter_rows()}
    stats["random_n_proteins"] = int(len(set(fr["gene_a"]) | set(fr["gene_b"])))
    stats["random_len_ratio_median"] = float(np.median(L[hi_i] / L[lo_i]))

    # ---------------- floor_ident ----------------
    g = (e.sort(["md5_v4", "acc"]).group_by("md5_v4", maintain_order=True)
         .agg(pl.col("acc"), pl.col("taxid"), pl.col("len").first().alias("len"))
         .filter(pl.col("acc").list.len() > 1))
    stats["ident_md5_groups"] = g.height
    stats["ident_possible_pairs"] = int((g["acc"].list.len() * (g["acc"].list.len() - 1) // 2).sum())
    out = []
    for b, (lo, hi) in enumerate(BINS):
        gb = g.filter((pl.col("len") >= lo) & (pl.col("len") < hi))
        idx = rng_i.permutation(gb.height)[: a.n_ident_per_bin]
        for i in idx:
            accs = gb["acc"][int(i)].to_list(); taxs = gb["taxid"][int(i)].to_list()
            u, v = rng_i.choice(len(accs), 2, replace=False)
            u, v = sorted((int(u), int(v)), key=lambda k: accs[k])
            out.append(dict(gene_a=accs[u], gene_b=accs[v], pair_set="floor_ident", hog="NA", split="NA",
                            species_a=taxs[u], species_b=taxs[v], relation="identical", len_a_fc=gb["len"][int(i)],
                            len_b_fc=gb["len"][int(i)], bin=BIN_NAMES[b]))
        stats.setdefault("ident_pairs_by_bin", {})[BIN_NAMES[b]] = [gb.height, min(len(idx), gb.height)]
    fi = pl.DataFrame(out)
    stats["ident_same_species_frac"] = float((fi["species_a"] == fi["species_b"]).mean())
    fi.write_csv(os.path.join(F, "floor_ident_pairs.tsv"), separator="\t")
    # accession list for fetch
    allacc = pl.concat([fr.select(pl.col("gene_a").alias("acc")), fr.select(pl.col("gene_b").alias("acc")),
                        fi.select(pl.col("gene_a").alias("acc")), fi.select(pl.col("gene_b").alias("acc"))]).unique()
    lk = pl.read_parquet(os.path.join(ROOT, "data/interim/panel/afdb_v4_lookup.parquet")).filter(pl.col("frag") == 1)
    al = allacc.join(lk.select("acc", "fc_key", "entry"), on="acc", how="left")
    assert al["fc_key"].null_count() == 0
    al.write_csv(os.path.join(F, "accessions.tsv"), separator="\t")
    stats["n_proteins_to_fetch"] = al.height
    json.dump(stats, open(os.path.join(F, "pair_build_stats.json"), "w"), indent=1)
    print(json.dumps(stats, indent=1))


if __name__ == "__main__":
    main()
