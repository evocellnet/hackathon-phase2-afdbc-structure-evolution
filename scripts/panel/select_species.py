#!/usr/bin/env python3
"""Select the species panel (UniProt reference proteomes) for the sequence-vs-structure divergence study.

Two independent panels (PI decision, design review round 1):
  --set euk   ~170 eukaryotes -> data/processed/panel/species_panel.tsv
  --set prok  ~100 prokaryotes (60 bacteria, 40 archaea) -> data/processed/panel/species_panel_prok.tsv

Eukaryote algorithm (deterministic, seed 20261030)
  1. Pool = UniProt reference proteomes (current release; metadata from the proteomes REST API, saved in
     data/raw/proteomes/meta/refprot_meta.tsv).  One proteome per NCBI species (QfO proteome preferred,
     then best quality).
  2. Filters
     * BUSCO complete >= 80% for Metazoa, Fungi, Viridiplantae; relaxed to >= 50% for the under-sampled
       protist supergroups (whose BUSCO eukaryota_odb scores are systematically low), and >= 15% for
       Metamonada (highly reduced parasite genomes, e.g. Giardia 24%).
     * AFDB v4 coverage: fraction of canonical accessions with an AFDB v4 model whose sequence equals the
       current UniProt sequence (exact md5 match; computed by map_proteomes_to_afdb.py once a proteome has
       been downloaded) >= --min-exact (default 0.90).  Before download a proxy is used: number of AFDB v4
       entries with the proteome's taxId / canonical protein count >= --min-exact.  The selection is
       iterated (select -> download -> map -> select) until stable.
     * Divergence time available in TimeTree 6 (pairwise API).
  3. Mandatory: QfO-2022 eukaryotes (Ancestra guide tree) and close sister-species pairs from SISTER_PAIRS
     (when they pass the filters).
  4. Breadth: one species per major clade (NCBI phylum; class within Chordata, Arthropoda, Ascomycota,
     Basidiomycota, Streptophyta, Chlorophyta, Apicomplexa) per supergroup (within its quota).
  5. Greedy fill up to --n-euk (at most 3 species per genus, sister pairs excepted).  Adding species c to a set S creates one new divergence node whose age is
     a(c) = min_{s in S} age(c, s).  We pick, at every step, the candidate whose new node age falls into the
     log10-age bin (1 Mya .. 2 Gya, width 0.3) that is most under-filled relative to a uniform target,
     combined with a supergroup quota term, a small quality term (BUSCO, coverage) and a seeded random
     tie-breaker.  age(c, s): TimeTree pairwise age for congeneric pairs (queried exactly; <= 12 species per genus);
     otherwise the age of their NCBI MRCA node (median TimeTree age of up to 3 representative pairs across
     the node's child clades).

Prokaryote algorithm: QfO bacteria/archaea + case-study species (kinases, actin/tubulin homologues) +
  phylum-quota round robin by quality (BUSCO, AFDB v4 exact coverage); seed 20261030.
"""
from __future__ import annotations

import argparse
import math
import os
import random
import re
import sys
from collections import Counter, defaultdict

import polars as pl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import ROOT, Taxonomy, TimeTree, parse_busco  # noqa: E402

SEED = 20261030
META = os.path.join(ROOT, "data/raw/proteomes/meta/refprot_meta.tsv")
QFO_TREE = os.path.join(ROOT, "tools/ancestra-mvp-ai/inputs/qfo2022/guide_tree.nwk")
COUNTS = os.path.join(ROOT, "data/interim/panel/afdb_v4_taxid_counts.parquet")
COVERAGE = os.path.join(ROOT, "data/interim/panel/proteome_v4_coverage.tsv")
OUTDIR = os.path.join(ROOT, "data/processed/panel")

# Close (and intermediate) sister-species pairs, requested explicitly (included when available + passing).
SISTER_PAIRS = [
    ("Homo sapiens", "Pan troglodytes"), ("Pan troglodytes", "Pan paniscus"), ("Homo sapiens", "Gorilla gorilla"),
    ("Mus musculus", "Rattus norvegicus"), ("Mus musculus", "Mus caroli"), ("Mus musculus", "Mus spretus"),
    ("Macaca mulatta", "Macaca fascicularis"),
    ("Drosophila melanogaster", "Drosophila simulans"), ("Drosophila simulans", "Drosophila sechellia"),
    ("Drosophila melanogaster", "Drosophila yakuba"), ("Drosophila pseudoobscura", "Drosophila persimilis"),
    ("Anopheles gambiae", "Anopheles arabiensis"), ("Anopheles gambiae", "Anopheles coluzzii"),
    ("Saccharomyces cerevisiae", "Saccharomyces paradoxus"), ("Saccharomyces cerevisiae", "Saccharomyces mikatae"),
    ("Candida albicans", "Candida dubliniensis"), ("Neurospora crassa", "Neurospora tetrasperma"),
    ("Aspergillus fumigatus", "Aspergillus fischeri"), ("Schizosaccharomyces pombe", "Schizosaccharomyces japonicus"),
    ("Cryptococcus neoformans", "Cryptococcus gattii"),
    ("Arabidopsis thaliana", "Arabidopsis lyrata"), ("Brassica rapa", "Brassica oleracea"),
    ("Oryza sativa", "Oryza glaberrima"), ("Oryza sativa", "Oryza brachyantha"), ("Solanum lycopersicum", "Solanum pennellii"),
    ("Caenorhabditis elegans", "Caenorhabditis briggsae"), ("Caenorhabditis briggsae", "Caenorhabditis nigoni"),
    ("Caenorhabditis remanei", "Caenorhabditis latens"), ("Brugia malayi", "Brugia pahangi"),
    ("Schistosoma mansoni", "Schistosoma haematobium"),
    ("Plasmodium falciparum", "Plasmodium reichenowi"), ("Plasmodium vivax", "Plasmodium knowlesi"),
    ("Plasmodium berghei", "Plasmodium yoelii"), ("Toxoplasma gondii", "Hammondia hammondi"),
    ("Toxoplasma gondii", "Neospora caninum"), ("Cryptosporidium parvum", "Cryptosporidium hominis"),
    ("Theileria annulata", "Theileria parva"),
    ("Leishmania major", "Leishmania infantum"), ("Leishmania infantum", "Leishmania donovani"),
    ("Trypanosoma brucei", "Trypanosoma cruzi"), ("Naegleria gruberi", "Naegleria fowleri"),
    ("Phytophthora infestans", "Phytophthora sojae"), ("Phytophthora sojae", "Phytophthora ramorum"),
    ("Entamoeba histolytica", "Entamoeba dispar"), ("Dictyostelium discoideum", "Dictyostelium purpureum"),
    ("Chlamydomonas reinhardtii", "Volvox carteri"), ("Ostreococcus tauri", "Ostreococcus lucimarinus"),
    ("Micromonas pusilla", "Micromonas commoda"), ("Giardia intestinalis", "Giardia muris"),
    ("Thalassiosira pseudonana", "Thalassiosira oceanica"), ("Galdieria sulphuraria", "Galdieria partita"),
]

# Eukaryote supergroup assignment (first lineage name that matches wins) and quotas (max species).
EUK_GROUPS = [  # (group, lineage token, quota, busco_min)
    ("Metazoa", "Metazoa", 66, 80), ("Fungi", "Fungi", 40, 80),
    ("Holozoa_other", "Opisthokonta", 4, 50),  # choanoflagellates, ichthyosporeans, filastereans, nucleariids
    ("Amoebozoa", "Amoebozoa", 6, 50), ("Apusozoa_CRuMs", "Apusozoa", 1, 50), ("CRuMs", "CRuMs", 1, 50),
    ("Land_plants", "Embryophyta", 20, 80), ("Charophytes", "Streptophyta", 2, 80),
    ("Chlorophyta", "Chlorophyta", 8, 80), ("Rhodophyta", "Rhodophyta", 4, 50),
    ("Glaucophyta", "Glaucocystophyceae", 1, 50),
    ("Stramenopiles", "Stramenopiles", 8, 50), ("Alveolata", "Alveolata", 10, 50), ("Rhizaria", "Rhizaria", 2, 50),
    ("Discoba", "Discoba", 8, 50), ("Metamonada", "Metamonada", 3, 15),
    ("Haptista", "Haptista", 2, 50), ("Cryptista", "Cryptophyceae", 1, 50),
]
# Documented BUSCO overrides (UniProt metadata artefacts; the genome is well annotated)
BUSCO_OVERRIDE = {"UP000008549": "C. briggsae AF16: UniProt reports BUSCO C:0.9%, an artefact (21.7k genes, "
                                 "96% exact AFDB v4 coverage, Ensembl/WormBase reference genome)"}
GENUS_CAP = 12
GENUS_MAX = 3        # greedy step: at most 3 species per genus (sister pairs may exceed it)
BREADTH_FRAC = 1.0   # breadth step may fill a supergroup quota completely
SPLIT_PHYLA = {"Chordata", "Arthropoda", "Ascomycota", "Basidiomycota", "Streptophyta", "Chlorophyta", "Apicomplexa"}
BIN_LO, BIN_HI, BIN_W = 0.0, math.log10(2000.0), 0.3


def load_pool(tx: Taxonomy, domain: str) -> pl.DataFrame:
    m = pl.read_csv(META, separator="\t", quote_char=None, infer_schema_length=0).rename(
        {"Proteome Id": "upid", "Organism": "name", "Organism Id": "taxid", "Protein count": "n_proteins",
         "BUSCO": "busco_str", "Taxonomic lineage": "lineage", "CPD": "cpd"})
    m = m.with_columns(pl.col("taxid").cast(pl.Int64), pl.col("n_proteins").cast(pl.Int64),
                       pl.col("busco_str").map_elements(parse_busco, return_dtype=pl.Float64).alias("busco"),
                       pl.col("lineage").str.split(", ").list.get(1, null_on_oob=True).alias("superkingdom"))
    m = m.filter(pl.col("superkingdom") == domain) if domain == "Eukaryota" else m.filter(
        pl.col("superkingdom").is_in(["Bacteria", "Archaea"]))
    cnt = pl.read_parquet(COUNTS)
    m = m.join(cnt, on="taxid", how="left").with_columns(pl.col("n_afdb_v4").fill_null(0))
    m = m.with_columns((pl.col("n_afdb_v4") / pl.col("n_proteins")).clip(upper_bound=1.0).alias("cov_proxy"))
    if os.path.exists(COVERAGE):
        cov = pl.read_csv(COVERAGE, separator="\t").select("upid", "n_canonical", "n_v4", "n_v4_exact")
        m = m.join(cov, on="upid", how="left")
        m = m.with_columns((pl.col("n_v4_exact") / pl.col("n_canonical")).alias("cov_exact"))
    else:
        m = m.with_columns(pl.lit(None, pl.Int64).alias(c) for c in ("n_canonical", "n_v4", "n_v4_exact")).with_columns(
            pl.lit(None, pl.Float64).alias("cov_exact"))
    m = m.with_columns(pl.coalesce("cov_exact", "cov_proxy").alias("cov"),
                       pl.col("cov_exact").is_not_null().alias("cov_is_exact"))
    m = m.with_columns(pl.col("taxid").map_elements(tx.species_of, return_dtype=pl.Int64).alias("species_taxid"))
    return m


def qfo_upids() -> set[str]:
    return {x.split("_")[0] for x in re.findall(r"UP\d+_\d+", open(QFO_TREE).read())}


def euk_group(lineage: str) -> tuple[str, int, int] | tuple[None, None, None]:
    toks = set(lineage.split(", "))
    for g, tok, q, b in EUK_GROUPS:
        if tok in toks:
            return g, q, b
    return None, None, None


def age_of(rec: dict | None) -> float | None:
    if not rec:
        return None
    a = rec.get("adjusted_age")
    if a is not None and a > 0:
        return a
    a = rec.get("precomputed_age")
    return a if a is not None and a > 0 else None


def write_euk_candidates(a, tx: Taxonomy):
    """Proteomes to download before selection so that the AFDB v4 exact-match coverage is known for every
    plausible candidate: BUSCO-passing (or QfO) proteomes with taxId-proxy coverage >= 0.5, plus every
    proteome of a species named in SISTER_PAIRS (the taxId proxy misses proteomes filed under another taxId)."""
    m = load_pool(tx, "Eukaryota")
    qfo = qfo_upids()
    names = {n for p in SISTER_PAIRS for n in p}
    keep = []
    for r in m.iter_rows(named=True):
        g, q, bmin = euk_group(r["lineage"])
        if g is None:
            continue
        spn = tx.name.get(r["species_taxid"], "")
        busco_ok = r["busco"] is not None and r["busco"] >= bmin
        if r["upid"] in qfo or spn in names or (busco_ok and r["cov_proxy"] >= 0.5):
            keep.append({"upid": r["upid"], "taxid": r["taxid"], "superkingdom": r["superkingdom"], "name": r["name"]})
    df = pl.DataFrame(keep)
    df.write_csv(os.path.join(ROOT, "data/interim/panel/euk_candidates.tsv"), separator="\t")
    print(f"[euk] {df.height} candidate proteomes to download")


def select_euk(a, tx: Taxonomy, tt: TimeTree):
    rng = random.Random(SEED)
    m = load_pool(tx, "Eukaryota")
    qfo = qfo_upids()
    rows = []
    for r in m.iter_rows(named=True):
        g, q, bmin = euk_group(r["lineage"])
        r["group"], r["busco_min"] = g, bmin
        r["is_qfo"] = r["upid"] in qfo
        rows.append(r)
    # filters
    stats = Counter()
    pool = []
    for r in rows:
        if r["group"] is None:
            stats["no_group"] += 1
            continue
        if not r["is_qfo"] and r["upid"] not in BUSCO_OVERRIDE and (r["busco"] is None or r["busco"] < r["busco_min"]):
            stats["busco"] += 1
            continue
        if r["cov"] is None or r["cov"] < a.min_exact:
            stats["afdb_cov"] += 1
            continue
        pool.append(r)
    # one proteome per species
    best = {}
    for r in sorted(pool, key=lambda r: (not r["is_qfo"], -(r["busco"] or 0) * r["cov"], r["upid"])):
        best.setdefault(r["species_taxid"], r)
    pool = sorted(best.values(), key=lambda r: r["upid"])
    print(f"[euk] metadata rows {len(rows)}; filtered out {dict(stats)}; pool (1/species) {len(pool)}", flush=True)
    sp = [r["species_taxid"] for r in pool]
    lin = {t: tx.lineage(t) for t in sp}

    # --- NCBI MRCA nodes of the pool and their TimeTree ages
    def mrca(x, y):
        la, lb = lin[x], lin[y]
        k = 0
        while k < min(len(la), len(lb)) and la[k] == lb[k]:
            k += 1
        return la[k - 1]
    children = defaultdict(lambda: defaultdict(list))  # node -> child -> [species]
    for t in sp:
        L = lin[t]
        for i in range(len(L) - 1):
            children[L[i]][L[i + 1]].append(t)
    nodes = {n: ch for n, ch in children.items() if len(ch) >= 2}
    q = []
    node_pairs = {}
    for n, ch in sorted(nodes.items()):
        keys = sorted(ch)
        pairs = []
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                pairs.append((min(ch[keys[i]]), min(ch[keys[j]])))
        rng.shuffle(pairs)
        node_pairs[n] = pairs[:3]
        q += pairs[:3]
    # exact pairwise ages for congeneric pairs (NCBI genus nodes are often polytomies with very different
    # pair ages); at most GENUS_CAP species per genus (seeded sample, QfO/sister-pair species first) are used.
    genus = {t: tx.rank_ancestor(t, "genus") for t in sp}
    priority_names = {n for p in SISTER_PAIRS for n in p}
    bygen = defaultdict(list)
    for t in sp:
        if genus[t]:
            bygen[genus[t]].append(t)
    for g_, ts in bygen.items():
        ts.sort(key=lambda t: (tx.name.get(t, "") not in priority_names, rng.random()))
        del ts[GENUS_CAP:]
    same_fam = [(x, y) for ts in bygen.values() for i, x in enumerate(ts) for y in ts[i + 1:]]
    print(f"[euk] {len(nodes)} MRCA nodes, {len(q)} node queries, {len(same_fam)} same-family pair queries", flush=True)
    res = tt.pairwise_many(q + same_fam)
    node_age = {}
    for n, prs in node_pairs.items():
        ages = sorted(x for x in (age_of(res.get((min(p), max(p)))) for p in prs) if x is not None)
        if ages:
            node_age[n] = ages[len(ages) // 2]

    def age(x, y):
        v = age_of(res.get((min(x, y), max(x, y))))
        if v is not None:
            return v
        n = mrca(x, y)
        while n not in node_age and n != 1:   # fall back to the closest dated ancestor
            n = tx.parent.get(n, 1)
        return node_age.get(n)

    # TimeTree availability: every species is queried against its closest pool member (and, if that fails,
    # against a second, more distant member); a species is kept if TimeTree dates either pair.
    dated = set()
    for (x, y), v in res.items():
        if age_of(v) is not None:
            dated.update((x, y))
    for attempt in range(2):
        chk = []
        for t in sp:
            if t in dated:
                continue
            L = lin[t]
            partners = []
            for depth in range(len(L) - 2, -1, -1):
                ch = children.get(L[depth], {})
                others = sorted(x for c, xs in ch.items() if c != L[depth + 1] for x in xs)
                if others:
                    partners.append(others[0])
                if len(partners) > attempt:
                    break
            if len(partners) > attempt:
                chk.append((t, partners[attempt]))
        r2 = tt.pairwise_many(chk)
        for (x, y) in chk:
            if age_of(r2.get((min(x, y), max(x, y)))) is not None:
                dated.add(x)
        res.update(r2)
    pool = [r for r in pool if r["species_taxid"] in dated]
    print(f"[euk] pool with TimeTree answer: {len(pool)}", flush=True)
    by = {r["species_taxid"]: r for r in pool}

    # --- mandatory species
    name2sp = {}
    for r in pool:
        nm = tx.name.get(r["species_taxid"], "")
        name2sp[nm] = r["species_taxid"]
    sel, reason = [], {}

    def add(t, why):
        if t not in reason:
            sel.append(t)
            reason[t] = why
        else:
            reason[t] += ";" + why
    for r in sorted(pool, key=lambda r: r["upid"]):
        if r["is_qfo"]:
            add(r["species_taxid"], "QfO2022")
    missing_pairs = []
    for x, y in SISTER_PAIRS:
        if x in name2sp and y in name2sp:
            add(name2sp[x], f"sister_pair:{x}|{y}")
            add(name2sp[y], f"sister_pair:{x}|{y}")
        else:
            missing_pairs.append((x, y, x in name2sp, y in name2sp))
    # --- breadth: one species per major clade (phylum; class within large phyla) per supergroup, up to
    # BREADTH_FRAC of the group quota, best quality first (seeded jitter breaks ties)
    quota = {g: qq for g, _, qq, _ in EUK_GROUPS}

    def clade_key(t):
        ph = tx.rank_ancestor(t, "phylum")
        if ph is None or tx.name.get(ph) in SPLIT_PHYLA:
            return tx.rank_ancestor(t, "class") or ph or tx.rank_ancestor(t, "order")
        return ph
    jit0 = {r["species_taxid"]: rng.random() for r in pool}
    gnow = Counter(by[t]["group"] for t in sel)
    have_keys = {clade_key(t) for t in sel}
    for g, _, qq, _ in EUK_GROUPS:
        members = [r for r in pool if r["group"] == g and r["species_taxid"] not in reason]
        keys = Counter(clade_key(r["species_taxid"]) for r in members)
        for k, _n in sorted(keys.items(), key=lambda kv: (-kv[1], kv[0] or 0)):
            if k in have_keys or gnow[g] >= max(1, int(BREADTH_FRAC * qq)):
                continue
            cands = sorted((r for r in members if clade_key(r["species_taxid"]) == k),
                           key=lambda r: (-(r["busco"] or 0) / 100.0 * r["cov"] - 1e-3 * jit0[r["species_taxid"]]))
            add(cands[0]["species_taxid"], f"breadth:{tx.name.get(k, k)}")
            have_keys.add(k)
            gnow[g] += 1
    # --- greedy log-uniform fill
    nb = math.ceil((BIN_HI - BIN_LO) / BIN_W)

    def binof(v):
        return min(nb - 1, max(0, int((math.log10(max(v, 1.0)) - BIN_LO) / BIN_W)))
    counts = [0] * nb
    order = list(sel)
    for i, t in enumerate(order[1:], 1):
        v = min((age(t, s) for s in order[:i] if age(t, s) is not None), default=None)
        if v is not None:
            counts[binof(v)] += 1
    gcount = Counter(by[t]["group"] for t in sel)
    gen_of = {r["species_taxid"]: tx.rank_ancestor(r["species_taxid"], "genus") for r in pool}
    gencount = Counter(gen_of[t] for t in sel)
    cand = [r["species_taxid"] for r in pool if r["species_taxid"] not in reason]
    best_age = {}
    for c in cand:
        vals = [v for v in (age(c, s) for s in sel) if v is not None]
        best_age[c] = min(vals) if vals else None
    jitter = {c: rng.random() for c in cand}
    while len(sel) < a.n_euk:
        target = (len(sel)) / nb
        best, bscore = None, -1e9
        for c in cand:
            if c in reason or best_age.get(c) is None:
                continue
            g = by[c]["group"]
            if gcount[g] >= quota[g] or (gen_of[c] and gencount[gen_of[c]] >= GENUS_MAX):
                continue
            b = binof(best_age[c])
            deficit = (target - counts[b]) / max(target, 1.0)
            need = (quota[g] - gcount[g]) / quota[g]
            qual = (by[c]["busco"] / 100.0) * by[c]["cov"]
            s = deficit + need + 0.1 * qual + 1e-3 * jitter[c]
            if s > bscore:
                best, bscore = c, s
        if best is None:
            break
        add(best, f"greedy:logbin{binof(best_age[best])}")
        counts[binof(best_age[best])] += 1
        gcount[by[best]["group"]] += 1
        gencount[gen_of[best]] += 1
        for c in cand:
            if c not in reason:
                v = age(c, best)
                if v is not None and (best_age[c] is None or v < best_age[c]):
                    best_age[c] = v
    print(f"[euk] selected {len(sel)}; group counts {dict(sorted(gcount.items()))}", flush=True)
    edges = [10 ** (BIN_LO + i * BIN_W) for i in range(nb + 1)]
    print("[euk] new-node ages per log bin (Mya):", ", ".join(f"{edges[i]:.0f}-{edges[i+1]:.0f}:{counts[i]}" for i in range(nb)))
    print(f"[euk] sister pairs not available/passing: {missing_pairs}")
    out = []
    for t in sel:
        r = by[t]
        out.append({"upid": r["upid"], "taxid": r["taxid"], "species_taxid": t, "name": r["name"],
                    "superkingdom": r["superkingdom"], "group": r["group"], "lineage": r["lineage"],
                    "busco": r["busco_str"], "busco_complete": r["busco"], "n_proteins": r["n_proteins"],
                    "afdb_v4_cov": round(r["cov"], 4), "afdb_v4_cov_is_exact": r["cov_is_exact"],
                    "is_qfo": r["is_qfo"], "nearest_age_at_selection_mya": None, "reason_included": reason[t]})
    # nearest panel neighbour age (approximate, for the table)
    for o in out:
        o["nearest_age_at_selection_mya"] = round(min((v for v in (age(o["species_taxid"], s) for s in sel if s != o["species_taxid"]) if v is not None), default=float("nan")), 2)
    return pl.DataFrame(out), missing_pairs


# ----------------------------------------------------------------------------- prokaryotes
PROK_CASE = [  # case-study species for kinases and cytoskeletal proteins (name prefix match on proteome name)
    "Escherichia coli (strain K12)", "Bacillus subtilis (strain 168)", "Caulobacter vibrioides", "Myxococcus xanthus",
    "Streptomyces coelicolor", "Mycobacterium tuberculosis (strain ATCC 25618", "Synechocystis sp. (strain ATCC 27184",
    "Nostoc sp. (strain PCC 7120", "Chlamydia trachomatis", "Gemmata obscuriglobus", "Prosthecobacter",
    "Thermus thermophilus", "Deinococcus radiodurans", "Aquifex aeolicus", "Thermotoga maritima",
    "Bacteroides thetaiotaomicron", "Rickettsia prowazekii", "Helicobacter pylori", "Staphylococcus aureus",
    "Haloferax volcanii", "Pyrobaculum calidifontis", "Thermoproteus tenax", "Sulfolobus acidocaldarius",
    "Saccharolobus solfataricus", "Nitrosopumilus maritimus", "Thermoplasma acidophilum", "Methanosarcina acetivorans",
    "Methanocaldococcus jannaschii", "Pyrococcus furiosus", "Halobacterium salinarum", "Nanoarchaeum equitans",
    "Candidatus Prometheoarchaeum syntrophicum", "Candidatus Korarchaeum cryptofilum", "Aeropyrum pernix",
    "Archaeoglobus fulgidus", "Methanobrevibacter smithii", "Thermococcus kodakarensis",
]


def select_prok(a, tx: Taxonomy):
    rng = random.Random(SEED)
    m = load_pool(tx, "prok")
    qfo = qfo_upids()
    out = []
    for dom, n_target, bmin in (("Bacteria", a.n_bact, 80.0), ("Archaea", a.n_arch, 50.0)):
        d = [r for r in m.filter(pl.col("superkingdom") == dom).iter_rows(named=True)]
        for r in d:
            L = r["lineage"].split(", ")
            # phylum-level label from NCBI ranks
            ph = tx.rank_ancestor(r["taxid"], "phylum")
            r["group"] = tx.name.get(ph, L[3] if len(L) > 3 else "unclassified") if ph else (L[3] if len(L) > 3 else "unclassified")
            r["is_qfo"] = r["upid"] in qfo
            asg = "Promethearchaeati" in L or "Asgard group" in L or "Promethearchaeota" in L
            dpann = "Nanobdellati" in L or "DPANN group" in L
            cpr = "Patescibacteriota" in L or "Bacteria candidate phyla" in L  # CPR: reduced genomes
            r["relaxed"] = asg or dpann or cpr
            r["min_cov"] = 0.5 if asg else a.min_exact  # Asgard: very few reference proteomes in 2021
        ok = [r for r in d if r["busco"] is not None and (r["busco"] >= (20.0 if r["relaxed"] else bmin))
              and r["cov"] is not None and r["cov"] >= r["min_cov"]]
        best = {}
        for r in sorted(ok, key=lambda r: (not r["is_qfo"], -(r["busco"] or 0) * r["cov"], r["upid"])):
            best.setdefault(r["species_taxid"], r)
        ok = sorted(best.values(), key=lambda r: r["upid"])
        print(f"[{dom}] metadata {len(d)}; passing {len(ok)}", flush=True)
        sel, reason = [], {}

        def add(r, why):
            if r["upid"] not in reason:
                sel.append(r)
                reason[r["upid"]] = why
        for r in ok:
            if r["is_qfo"]:
                add(r, "QfO2022")
        for pat in PROK_CASE:
            hits = sorted([r for r in ok if r["name"].startswith(pat)], key=lambda r: (-(r["busco"] or 0) * r["cov"], r["upid"]))
            if hits:
                add(hits[0], f"case_study:{pat}")
        # phylum round robin: phyla ordered by size of pool (desc), best quality first within phylum (+jitter)
        byph = defaultdict(list)
        for r in ok:
            if r["upid"] not in reason:
                byph[r["group"]].append(r)
        for ph in byph:
            byph[ph].sort(key=lambda r: (-(r["busco"] or 0) * r["cov"] - 1e-3 * rng.random()))
        have = Counter(r["group"] for r in sel)
        phyla = sorted(byph, key=lambda p: (have[p], -len(byph[p]), p))
        while len(sel) < n_target and any(byph.values()):
            phyla = sorted([p for p in byph if byph[p]], key=lambda p: (have[p], -len(byph[p]), p))
            p = phyla[0]
            # within a phylum, prefer a genus not yet selected
            gen_have = {tx.rank_ancestor(r["taxid"], "genus") for r in sel}
            pick = next((r for r in byph[p] if tx.rank_ancestor(r["taxid"], "genus") not in gen_have), byph[p][0])
            byph[p].remove(pick)
            add(pick, f"phylum_round_robin:{p}")
            have[p] += 1
        print(f"[{dom}] selected {len(sel)}: {dict(Counter(r['group'] for r in sel))}", flush=True)
        for r in sel:
            out.append({"upid": r["upid"], "taxid": r["taxid"], "species_taxid": r["species_taxid"], "name": r["name"],
                        "superkingdom": dom, "group": r["group"], "lineage": r["lineage"], "busco": r["busco_str"],
                        "busco_complete": r["busco"], "n_proteins": r["n_proteins"], "afdb_v4_cov": round(r["cov"], 4),
                        "afdb_v4_cov_is_exact": r["cov_is_exact"], "is_qfo": r["is_qfo"],
                        "reason_included": reason[r["upid"]]})
    return pl.DataFrame(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--set", choices=["euk", "prok", "euk-candidates"], required=True)
    ap.add_argument("--n-euk", type=int, default=170)
    ap.add_argument("--n-bact", type=int, default=60)
    ap.add_argument("--n-arch", type=int, default=40)
    ap.add_argument("--min-exact", type=float, default=0.90)
    ap.add_argument("--out")
    a = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)
    tx = Taxonomy()
    if a.set == "euk-candidates":
        write_euk_candidates(a, tx)
        return
    if a.set == "euk":
        df, _ = select_euk(a, tx, TimeTree())
        out = a.out or os.path.join(OUTDIR, "species_panel.tsv")
    else:
        df = select_prok(a, tx)
        out = a.out or os.path.join(OUTDIR, "species_panel_prok.tsv")
    df = df.with_columns(pl.col("upid").alias("leaf"))
    df = df.with_columns((pl.col("upid") + "_" + pl.col("taxid").cast(pl.Utf8)).alias("leaf"))
    df.write_csv(out, separator="\t")
    print(f"wrote {out} ({df.height} rows)")


if __name__ == "__main__":
    main()
