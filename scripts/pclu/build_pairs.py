#!/usr/bin/env python3
"""P-clu stage 1: build one analysis-ready parquet of member->representative alignments.

Inputs (all under data/raw/):
  afdb_clusters/v3/8-cluster_alignment/aln-*.tsv   (no header; columns verified: query=representative,
        target=member, fident alnlen mismatch gapopen qstart qend tstart tend evalue bits lddt alntmscore)
  afdb_clusters/v3/1-AFDBClusters-entryId_repId_taxId.tsv.gz  (memberId, repId, taxId)
  afdb_clusters/v3/2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv.gz
  ncbi_taxonomy/taxdump/{nodes,merged}.dmp, ncbi_taxonomy/new_taxdump/rankedlineage.dmp
Output: data/processed/pclu/pairs.parquet  (+ results/pclu/build_checks.json)
"""
import json, sys, time
from pathlib import Path
import numpy as np
import polars as pl

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data/raw"
V3 = RAW / "afdb_clusters/v3"
OUT = ROOT / "data/processed/pclu"
RES = ROOT / "results/pclu"
OUT.mkdir(parents=True, exist_ok=True); RES.mkdir(parents=True, exist_ok=True)
ALN = next((V3 / "8-cluster_alignment").glob("aln-*.tsv"))

STD_RANKS = ["species", "genus", "family", "order", "class", "phylum", "kingdom", "domain"]
LIN_RANKS = STD_RANKS  # rankedlineage columns
DOMAIN_IDS = {2: "Bacteria", 2157: "Archaea", 2759: "Eukaryota", 10239: "Viruses"}


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def read_dmp(path, ncols):
    rows = []
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t|\t")
            f[-1] = f[-1].rstrip("\t|")
            rows.append(f[:ncols])
    return rows


def load_taxonomy():
    """Returns dict of numpy arrays indexed by taxid (size max+1): parent, rank_code, domain_code,
    std_rank_of_self_or_lowest_ancestor, plus merged map and lineage name codes."""
    log("taxonomy: nodes/merged")
    nodes = read_dmp(RAW / "ncbi_taxonomy/taxdump/nodes.dmp", 3)
    merged = read_dmp(RAW / "ncbi_taxonomy/taxdump/merged.dmp", 2)
    ids = np.array([int(r[0]) for r in nodes]); par = np.array([int(r[1]) for r in nodes])
    rank = [r[2] for r in nodes]
    N = max(ids.max(), max(int(m[0]) for m in merged)) + 1
    parent = np.zeros(N, dtype=np.int64); parent[ids] = par
    present = np.zeros(N, dtype=bool); present[ids] = True
    rank_arr = np.full(N, "", dtype=object); rank_arr[ids] = rank
    merged_to = np.arange(N, dtype=np.int64)
    for o, n in merged:
        merged_to[int(o)] = int(n)
    # domain code by walking up (vectorised, max depth ~ 60)
    dom = np.zeros(N, dtype=np.int8)  # 0 other/unknown, 1 Bacteria, 2 Archaea, 3 Eukaryota, 4 Viruses
    dmap = {2: 1, 2157: 2, 2759: 3, 10239: 4}
    cur = np.arange(N, dtype=np.int64)
    for _ in range(80):
        for tid, code in dmap.items():
            dom[(cur == tid) & (dom == 0)] = code
        nxt = parent[cur]
        nxt[~present[cur]] = 1
        cur = nxt
    # lowest standard-rank ancestor-or-self (index into STD_RANKS; 99 = none, i.e. root/cellular organisms)
    std_idx = np.full(N, 99, dtype=np.int8)
    r2i = {r: i for i, r in enumerate(STD_RANKS)}
    own = np.array([r2i.get(r, -1) if present[i] else -1 for i, r in enumerate(rank_arr)])
    # species-level or below: treat 'strain','subspecies','varietas','forma', 'isolate', 'serotype' etc as species
    below_species = {"subspecies", "strain", "varietas", "forma", "forma specialis", "isolate", "serotype",
                    "serogroup", "biotype", "pathogroup", "morph", "genotype"}
    for i, r in enumerate(rank_arr):
        if present[i] and r in below_species:
            own[i] = 0
    cur = np.arange(N, dtype=np.int64)
    for _ in range(80):
        m = (std_idx == 99) & (own[cur] >= 0)
        std_idx[m] = own[cur][m]
        nxt = parent[cur]; nxt[~present[cur]] = 1; cur = nxt
    log("taxonomy: rankedlineage")
    rl = read_dmp(RAW / "ncbi_taxonomy/new_taxdump/rankedlineage.dmp", 10)
    # columns: taxid, name, species, genus, family, order, class, phylum, kingdom, domain
    names = {}
    L = np.zeros((N, len(LIN_RANKS)), dtype=np.int32)  # 0 = empty
    for r in rl:
        t = int(r[0])
        for j in range(len(LIN_RANKS)):
            nm = r[2 + j].strip()
            if nm:
                L[t, j] = names.setdefault(nm, len(names) + 1)
    return dict(parent=parent, present=present, rank=rank_arr, merged_to=merged_to, dom=dom,
                std_idx=std_idx, L=L, N=N)


def main():
    t0 = time.time()
    tax = load_taxonomy()
    log("reading file 8 (alignments)")
    cols = ["query", "target", "fident", "alnlen", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend",
            "evalue", "bits", "lddt", "tm"]
    schema = {c: pl.Float64 for c in cols}
    schema.update(query=pl.Utf8, target=pl.Utf8, alnlen=pl.Int32, mismatch=pl.Int32, gapopen=pl.Int32,
                  qstart=pl.Int32, qend=pl.Int32, tstart=pl.Int32, tend=pl.Int32)
    df = pl.read_csv(ALN, separator="\t", has_header=False, new_columns=cols, schema_overrides=schema)
    strip = lambda c: pl.col(c).str.replace(r"^AF-", "").str.replace(r"-F1-model_v3\.cif$", "").alias(c)
    df = df.with_columns(strip("query"), strip("target"))
    log("rows", df.height)
    checks = {"n_rows_file8": df.height}
    checks["n_exact_duplicate_pairs"] = int(df.height - df.select("query", "target").unique().height)
    checks["n_nan"] = {c: int(df[c].is_null().sum() + df[c].is_nan().sum()) for c in ["fident", "evalue", "lddt", "tm"]}

    log("file 1 membership")
    mem = pl.read_csv(V3 / "1-AFDBClusters-entryId_repId_taxId.tsv.gz", separator="\t", has_header=False,
                      new_columns=["target", "query_f1", "taxid_member"],
                      schema_overrides={"target": pl.Utf8, "query_f1": pl.Utf8, "taxid_member": pl.Int64})
    clu = pl.read_csv(V3 / "2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv.gz", separator="\t",
                      has_header=False, new_columns=["query", "isDark", "nMem", "repLen", "avgLen", "repPlddt",
                                                     "avgPlddt", "lca"],
                      schema_overrides={"query": pl.Utf8})
    # verify orientation: file-8 (query,target) must equal file-1 (repId, memberId)
    j = df.join(mem, on="target", how="left")
    checks["orientation_query_is_rep_frac"] = float((j["query"] == j["query_f1"]).mean())
    checks["n_target_not_in_file1"] = int(j["query_f1"].is_null().sum())
    del j
    checks["n_query_not_in_file2"] = int(df.select("query").unique().join(clu.select("query"), on="query", how="anti").height)
    # rep taxid
    reptax = mem.filter(pl.col("target") == pl.col("query_f1")).select(pl.col("target").alias("query"),
                                                                       pl.col("taxid_member").alias("taxid_rep"))
    df = (df.join(mem.select("target", "taxid_member"), on="target", how="left")
            .join(reptax, on="query", how="left")
            .join(clu, on="query", how="left"))
    df = df.with_columns((pl.col("query") == pl.col("target")).alias("is_self"))
    checks["n_self_pairs"] = int(df["is_self"].sum())
    checks["n_clusters_file2"] = clu.height

    # taxonomy annotations (resolve merged ids)
    def tx(col):
        a = df[col].fill_null(0).to_numpy().astype(np.int64)
        a = np.where(a < tax["N"], a, 0)
        return tax["merged_to"][a]
    tm_ = tx("taxid_member"); tr_ = tx("taxid_rep"); lca_ = tx("lca")
    dnames = np.array(["Other", "Bacteria", "Archaea", "Eukaryota", "Viruses"])
    d_m = dnames[tax["dom"][tm_]]; d_r = dnames[tax["dom"][tr_]]
    # pair MRCA rank: lowest standard rank at which both lineages share the same non-empty name;
    # 0=same species name ... 7=same domain label ; 8 = share none (different domain) ; 9 = same taxid
    Lm = tax["L"][tm_]; Lr = tax["L"][tr_]
    eq = (Lm == Lr) & (Lm > 0)
    first = np.where(eq.any(axis=1), eq.argmax(axis=1), 8)
    first = np.where(tm_ == tr_, 9, first).astype(np.int8)
    # a species-level shared name in a higher rank implies lower; argmax finds the lowest shared rank
    pair_rank_names = ["species", "genus", "family", "order", "class", "phylum", "kingdom", "domain", "none", "same_taxid"]
    # fallback: if rankedlineage gave no shared name but the domain codes agree, call it 'domain'
    same_dom = (tax["dom"][tm_] == tax["dom"][tr_]) & (tax["dom"][tm_] > 0)
    first = np.where((first == 8) & same_dom, 7, first).astype(np.int8)
    lca_std = tax["std_idx"][lca_]
    lca_rank_name = np.array(STD_RANKS + ["root_or_cellular"], dtype=object)[np.where(lca_std == 99, 8, lca_std)]
    lca_rank_name = np.where(np.isin(lca_, [1, 131567]), "root_or_cellular", lca_rank_name)
    lca_dom = dnames[tax["dom"][lca_]]
    df = df.with_columns(
        pl.Series("taxid_member", tm_), pl.Series("taxid_rep", tr_), pl.Series("lca", lca_),
        pl.Series("dom_member", d_m), pl.Series("dom_rep", d_r),
        pl.Series("pair_mrca", np.array(pair_rank_names, dtype=object)[first]),
        pl.Series("lca_rank", lca_rank_name), pl.Series("lca_dom", lca_dom),
    )
    df = df.with_columns(
        (pl.col("dom_member") == pl.col("dom_rep")).alias("same_domain"),
        ((pl.col("qend") - pl.col("qstart") + 1) / pl.col("repLen")).alias("qcov"),
        pl.col("isDark").cast(pl.Boolean),
    )
    checks["frac_member_taxid_unresolved"] = float((tm_ == 0).mean())
    df.write_parquet(OUT / "pairs.parquet", compression="zstd")
    (RES / "build_checks.json").write_text(json.dumps(checks, indent=1))
    log("done", df.height, "rows;", f"{time.time()-t0:.0f}s", json.dumps(checks))


if __name__ == "__main__":
    main()
