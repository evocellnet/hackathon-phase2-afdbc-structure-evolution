#!/usr/bin/env python3
"""HOG construction (ANALYSIS_PLAN §1.6, §1.1; E6 triplets §4.6).

  python scripts/hogs/build_hogs.py --orthoxml data/interim/ancestra_euk180/ancestra/hogs.orthoxml --threads 48

Steps (every step logs counts in -> out to <logs>/hogs_*.tsv):
 1. Gene-model filters on the members of every top-level Ancestra HOG:
      md5(UniProt 2026) == md5(AFDB v4)  [protein_table.seq_match; genes outside the panel fail here]
      not a fragment: UniProt '(Fragment)' or AFDB Clusters v3 cluFlag 3  [data/interim/hogs/gene_flags.parquet]
      mean pLDDT >= 70
      length within 0.5-2x the median length of the HOG members that passed the three filters above
    1:1 set = species with exactly one filtered copy; keep HOGs with >= 8 (n_species_1to1).
 2. Presample <= 60 species by greedy panel-node coverage (gain = panel-tree ancestors not yet covered; ties by a
    seeded draw, stream hog_presample, key crc32(hog)).
 3. Clade-closed selection: panel internal nodes C ordered by (#presample species in C desc, older first, node id);
    S = presample ∩ C; accept the first C with |S| >= 8 and every pair aa_ml_seq <= 2 (NaN/saturated = fail);
    |S| > 30 -> 30 by greedy node coverage within C. aa_ml_seq is computed with the pairs pipeline's seqonly
    code path (scripts/hogs/seqdist.py; identical to the pipeline to 2e-15). Default = lazy evaluation (pairs are
    computed on demand, deepest species pairs first, and a clade is rejected at its first failing pair); the
    selection is identical to computing all presample pairs (--exhaustive), which is the plan's literal recipe.
 4. Sample <= 7,000 HOGs: equal quotas per root-age stratum (<100, 100-300, 300-1000, >=1000 Mya), unused quota
    redistributed proportionally to the other strata's spare HOGs (seeded, stream hog_sample); per stratum a seeded
    permutation, take a prefix. If primary + E6 pairs > budget (3.5M), all quotas are scaled down (prefixes, so
    the smaller sample is nested in the larger) until it fits.
 5. Split 20/80 within strata (seeded permutation, stream split) -> split.tsv, written ONCE; if it exists it is
    read and checked, never regenerated.
 6. Primary pairs: all pairs within S_h; labels from scripts/ancestra/hog_pair_labels.py (label_pairs).
 7. E6 (§4.6), from the filtered membership of sampled HOGs: P1,P2 = in_paralog pair, duplication level P not a
    leaf, aa_pid_seq <= 0.95 and aa_ml_seq <= 2; O1 = species of the sister clade(s) of P (guide tree of the
    OrthoXML) with exactly one filtered gene in the HOG, that gene ortholog to both copies; O2 = same for the
    sister(s) of parent(P); |O1| >= 1. P1 = seeded coin flip (stream e6_p1_choice).
    Pair sets: e6_para (P1,P2), e6_orth (P1|P2, O1), e6_ldo (P1|P2, O2).
Outputs: <outdir>/hogs.parquet, <pairlists>/{primary,e6}.tsv, <pairlists>/e6_triplets.tsv, split.tsv,
         <outdir>/hogs_build_summary.json
"""
from __future__ import annotations

import argparse
import itertools
import json
import multiprocessing as mp
import sys
import time
import zlib
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "tools/ancestra-mvp-ai"))
sys.path.insert(0, str(ROOT / "scripts/ancestra"))
import seqdist  # noqa: E402  (sets BLAS threads, imports strucevo)
import treeutil as T  # noqa: E402

SEED = 20261030
STAGES = {"hog_presample": 0, "hog_sample": 1, "split": 2, "e6_p1_choice": 14}
STRATA = [(1, 0, 100), (2, 100, 300), (3, 300, 1000), (4, 1000, np.inf)]


def rng_for(stage: str, *key) -> np.random.Generator:
    return np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(STAGES[stage],) + tuple(key)))


def crc(s: str) -> int:
    return zlib.crc32(s.encode())


def taxid_of(leaf: str) -> int:
    return int(leaf.rsplit("_", 1)[1])


class Log:
    def __init__(self, d: Path, tag: str):
        self.d = d; self.tag = tag; self.tables = defaultdict(list)
        d.mkdir(parents=True, exist_ok=True)

    def row(self, table, **kw):
        self.tables[table].append(kw)
        print(f"[{table}] " + " ".join(f"{k}={v}" for k, v in kw.items()), flush=True)

    def write(self):
        for t, rows in self.tables.items():
            pd.DataFrame(rows).to_csv(self.d / f"{self.tag}{t}.tsv", sep="\t", index=False)


# ------------------------------------------------------------------------------------------------ panel tree
class Panel:
    def __init__(self, nwk: Path, ages: Path):
        self.root = T.parse_newick(nwk.read_text())
        self.nodes, self.leafset = T.index_tree(self.root)
        na = pd.read_csv(ages, sep="\t")
        age = dict(zip(na.node, na["median"]))
        self.age = np.array([age.get(n.name, np.nan) for n in self.nodes])
        self.name = [n.name for n in self.nodes]
        self.leaf_id = {n.name: n.id for n in self.nodes if n.is_leaf()}
        self.anc = {}
        for n in self.nodes:
            if n.is_leaf():
                a, x = [], n
                while x is not None:
                    a.append(x.id); x = x.parent
                self.anc[n.name] = a
        self.internal = [n.id for n in self.nodes if not n.is_leaf()]

    def mrca(self, leaves) -> int:
        leaves = list(leaves)
        common = set(self.anc[leaves[0]])
        for l in leaves[1:]:
            common &= set(self.anc[l])
        return max(common, key=lambda i: -len(self.leafset[i]))  # smallest clade

    def pair_age(self, a, b) -> float:
        sa = set(self.anc[b])
        m = next(x for x in self.anc[a] if x in sa)
        v = self.age[m]
        return -1.0 if np.isnan(v) else v


def greedy_cover(panel: Panel, species: list, k: int, rng) -> list:
    """Greedy panel-node coverage: repeatedly add the species with the most not-yet-covered ancestor nodes."""
    if len(species) <= k:
        return sorted(species)
    rem = sorted(species)
    covered, out = set(), []
    while len(out) < k:
        gains = np.array([len(set(panel.anc[s]) - covered) for s in rem])
        best = np.flatnonzero(gains == gains.max())
        j = best[rng.integers(len(best))] if len(best) > 1 else best[0]
        s = rem.pop(j)
        out.append(s); covered |= set(panel.anc[s])
    return sorted(out)


class HogIndex:
    """Pair labels for one top-level HOG by the deepest-common-group rule, computed per pair on demand.
    Same semantics as scripts/ancestra/hog_pair_labels.py:label_pairs (checked in tests/test_hogs.py)."""

    def __init__(self, top, sp_of, par):
        from ancestra.orthoxml import Group
        self.sp_of, self.par = sp_of, par
        self.path = {}
        st = [(top, [top])]
        while st:
            g, p = st.pop()
            for c in g.children:
                if isinstance(c, Group):
                    st.append((c, p + [c]))
                else:
                    self.path[c] = p + [c]
        self._cnt = {}

    def cnt(self, node):
        k = id(node)
        if k not in self._cnt:
            if isinstance(node, str):
                self._cnt[k] = Counter([self.sp_of[node]])
            else:
                self._cnt[k] = Counter(self.sp_of[x] for x in node.iter_genes())
        return self._cnt[k]

    def label(self, x, y):
        """-> (relation, event_level, dup_branch_from, one_to_one, dcg_id)"""
        px, py = self.path[x], self.path[y]
        k = 0
        while k + 1 < min(len(px), len(py)) and px[k + 1] is py[k + 1]:
            k += 1
        g = px[k]; cx, cy = px[k + 1], py[k + 1]
        sa, sb = self.sp_of[x], self.sp_of[y]
        if g.kind == "ortholog":
            if sa == sb:
                return "inconsistent", g.level, "", 0, g.id or ""
            ci, cj = self.cnt(cx), self.cnt(cy)
            o2o = int(ci[sa] == 1 and cj[sb] == 1 and ci[sb] == 0 and cj[sa] == 0)
            return "ortholog", g.level, "", o2o, g.id or ""
        frm = self.par.get(g.level, "") if g.level else ""
        return ("in_paralog" if sa == sb else "out_paralog"), g.level or "", frm, 0, g.id or ""


# ------------------------------------------------------------------------------------------------ main
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--orthoxml", required=True)
    ap.add_argument("--outdir", default=str(ROOT / "data/interim"))
    ap.add_argument("--pairlists", default=None, help="default <outdir>/pairlists")
    ap.add_argument("--split", default=None, help="default <outdir>/split.tsv (written once, never regenerated)")
    ap.add_argument("--logs", default=str(ROOT / "logs/filters"))
    ap.add_argument("--log-tag", default="hogs_")
    ap.add_argument("--panel-tree", default=str(ROOT / "data/interim/timetree/panel_tree.nwk"))
    ap.add_argument("--node-ages", default=str(ROOT / "data/interim/timetree/node_ages.tsv"))
    ap.add_argument("--protein-table", default=str(ROOT / "data/processed/panel/protein_table.parquet"))
    ap.add_argument("--gene-flags", default=str(ROOT / "data/interim/hogs/gene_flags.parquet"))
    ap.add_argument("--fasta-dir", default=str(ROOT / "data/interim/ancestra_input"))
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--min-plddt", type=float, default=70.0)
    ap.add_argument("--min-1to1", type=int, default=8)
    ap.add_argument("--max-presample", type=int, default=60)
    ap.add_argument("--max-sel", type=int, default=30)
    ap.add_argument("--dmax", type=float, default=2.0)
    ap.add_argument("--alpha-aa", type=float, default=1.0, help="LG gamma alpha for aa_ml_seq (placeholder 1.0)")
    ap.add_argument("--n-hogs", type=int, default=7000)
    ap.add_argument("--pair-budget", type=float, default=3.5e6)
    ap.add_argument("--tune-frac", type=float, default=0.2)
    ap.add_argument("--pid-max", type=float, default=0.95)
    ap.add_argument("--exhaustive", action="store_true", help="compute aa_ml_seq for ALL presample pairs")
    ap.add_argument("--k0", type=int, default=48, help="lazy mode: first pair batch per HOG and clade")
    ap.add_argument("--e6-max-o", type=int, default=None,
                    help="cap |O1| and |O2| per triplet (seeded subset). Default: no cap (plan-literal)")
    ap.add_argument("--e6-max-per-dcg", type=int, default=None,
                    help="cap triplets per paralog group dcg_id (seeded subset). Default: no cap (plan-literal)")
    ap.add_argument("--no-cache", action="store_true", help="ignore/skip <outdir>/hogs/seqdist_cache.parquet")
    ap.add_argument("--projection-only", action="store_true",
                    help="stop after step 3 + E6 counting on all eligible HOGs; no split/pair lists written")
    a = ap.parse_args(argv)
    t0 = time.time()
    outdir = Path(a.outdir); outdir.mkdir(parents=True, exist_ok=True)
    pldir = Path(a.pairlists) if a.pairlists else outdir / "pairlists"
    split_path = Path(a.split) if a.split else outdir / "split.tsv"
    log = Log(Path(a.logs), a.log_tag)
    summary = {"args": vars(a)}

    # ---------------------------------------------------------------- load
    from ancestra.orthoxml import Group, read_orthoxml
    ox = read_orthoxml(a.orthoxml)
    tops = [g for g in ox.groups if isinstance(g, Group)]
    sp_of = {k: v.species for k, v in ox.genes.items()}
    print(f"orthoxml: {len(ox.genes)} genes, {len(tops)} top-level HOGs ({time.time()-t0:.0f}s)", flush=True)
    panel = Panel(Path(a.panel_tree), Path(a.node_ages))
    pt = pd.read_parquet(a.protein_table, columns=["acc", "leaf", "len_2026", "seq_match", "mean_plddt"])
    gf = pd.read_parquet(a.gene_flags)
    pt = pt.merge(gf[["acc", "frag_uniprot", "cluflag"]], on="acc", how="left")
    pt = pt.set_index("acc")
    hog_genes = {g.id: list(g.iter_genes()) for g in tops}
    all_genes = [x for v in hog_genes.values() for x in v]
    gi = pt.reindex(all_genes)
    in_panel = gi["leaf"].notna().to_numpy()
    f_md5 = in_panel & gi["seq_match"].fillna(False).to_numpy(bool)
    f_frag = ~(gi["frag_uniprot"].fillna(False).to_numpy(bool) | (gi["cluflag"].fillna(0).to_numpy() == 3))
    f_pl = gi["mean_plddt"].fillna(-1).to_numpy() >= a.min_plddt
    glen = gi["len_2026"].fillna(-1).to_numpy()
    gpl = dict(zip(all_genes, gi["mean_plddt"].to_numpy()))
    # per-gene pass flags in HOG order
    pos = 0
    passed = {}            # gene -> True if passes all filters
    hog_stats = {}
    cnt = Counter()
    for g in tops:
        gs = hog_genes[g.id]; n = len(gs); sl = slice(pos, pos + n); pos += n
        m1 = f_md5[sl]; m2 = m1 & f_frag[sl]; m3 = m2 & f_pl[sl]
        L = glen[sl]
        if m3.any():
            med = float(np.median(L[m3]))
            m4 = m3 & (L >= 0.5 * med) & (L <= 2.0 * med)
        else:
            med, m4 = np.nan, m3
        for x, ok in zip(gs, m4):
            if ok:
                passed[x] = True
        cnt["genes_in"] += n; cnt["md5"] += int(m1.sum()); cnt["frag"] += int(m2.sum())
        cnt["plddt"] += int(m3.sum()); cnt["length"] += int(m4.sum())
        cnt["not_in_panel"] += int((~in_panel[sl]).sum())
        for k, m in (("hog_md5", m1), ("hog_frag", m2), ("hog_plddt", m3), ("hog_length", m4)):
            cnt[k] += int(m.any())
        spc = Counter(sp_of[x] for x, ok in zip(gs, m4) if ok)
        one = sorted(s for s, c in spc.items() if c == 1)
        hog_stats[g.id] = dict(n_genes=n, n_genes_filtered=int(m4.sum()), median_len=med, one=one,
                               one_gene={s: next(x for x, ok in zip(gs, m4) if ok and sp_of[x] == s) for s in one})
    cand = [h for h, s in hog_stats.items() if len(s["one"]) >= a.min_1to1]
    nh = len(tops)
    log.row("genes", step="input", genes_in=cnt["genes_in"], genes_out=cnt["genes_in"], hogs_in=nh, hogs_out=nh)
    log.row("genes", step="md5_match", genes_in=cnt["genes_in"], genes_out=cnt["md5"], hogs_in=nh,
            hogs_out=cnt["hog_md5"], note=f"{cnt['not_in_panel']} genes not in panel protein table")
    log.row("genes", step="not_fragment", genes_in=cnt["md5"], genes_out=cnt["frag"], hogs_in=cnt["hog_md5"],
            hogs_out=cnt["hog_frag"])
    log.row("genes", step="mean_plddt_ge70", genes_in=cnt["frag"], genes_out=cnt["plddt"], hogs_in=cnt["hog_frag"],
            hogs_out=cnt["hog_plddt"])
    log.row("genes", step="length_0.5_2x_median", genes_in=cnt["plddt"], genes_out=cnt["length"],
            hogs_in=cnt["hog_plddt"], hogs_out=cnt["hog_length"])
    n1 = sum(len(s["one"]) for s in hog_stats.values())
    log.row("genes", step=f"1to1_species_ge{a.min_1to1}", genes_in=cnt["length"],
            genes_out=sum(len(hog_stats[h]["one"]) for h in cand), hogs_in=cnt["hog_length"], hogs_out=len(cand),
            note=f"{n1} 1:1 genes in all HOGs")
    summary["filters"] = dict(cnt)
    # panel species check
    unknown_sp = {sp_of[x] for x in passed} - set(panel.leaf_id)
    if unknown_sp:
        sys.exit(f"species not in panel tree: {sorted(unknown_sp)[:5]}")

    # ---------------------------------------------------------------- presample
    pres = {}
    for h in cand:
        pres[h] = greedy_cover(panel, hog_stats[h]["one"], a.max_presample, rng_for("hog_presample", crc(h)))
    log.row("selection", step="presample", hogs_in=len(cand), hogs_out=len(cand),
            note=f"mean presample size {np.mean([len(v) for v in pres.values()]):.1f}; "
                 f"total presample pairs {sum(len(v)*(len(v)-1)//2 for v in pres.values())}")

    # sequences of all filtered genes
    t1 = time.time()
    seqdist.SEQ.update(seqdist.read_fasta_dir(Path(a.fasta_dir), set(passed)))
    miss = set(passed) - set(seqdist.SEQ)
    if miss:
        sys.exit(f"{len(miss)} filtered genes without sequence, e.g. {sorted(miss)[:3]}")
    print(f"sequences: {len(seqdist.SEQ)} ({time.time()-t1:.0f}s)", flush=True)
    pool = mp.get_context("fork").Pool(a.threads)

    memo = {}             # frozenset{gene,gene} -> (ml, pid, sat)
    n_computed = [0]
    cache = outdir / "hogs" / f"seqdist_cache_alpha{a.alpha_aa:g}.parquet"
    n_cached = 0
    if cache.exists() and not a.no_cache:
        cdf = pd.read_parquet(cache)
        memo.update({frozenset((x, y)): (m, q, bool(z)) for x, y, m, q, z in
                     zip(cdf.a, cdf.b, cdf.ml, cdf.pid, cdf.sat)})
        n_cached = len(memo)
        print(f"loaded {n_cached} cached seqonly distances from {cache}", flush=True)

    def save_cache():
        if a.no_cache or len(memo) == n_cached:
            return
        cache.parent.mkdir(parents=True, exist_ok=True)
        ks = list(memo)
        pd.DataFrame({"a": [min(k) for k in ks], "b": [max(k) if len(k) > 1 else min(k) for k in ks],
                      "ml": [memo[k][0] for k in ks], "pid": [memo[k][1] for k in ks],
                      "sat": [memo[k][2] for k in ks]}).to_parquet(cache, index=False)

    def compute(pairs):
        pairs = [p for p in dict.fromkeys(tuple(sorted(p)) for p in pairs) if frozenset(p) not in memo]
        if not pairs:
            return
        ml, pid, sat = seqdist.seq_distances(pairs, a.threads, a.alpha_aa, pool=pool)
        for p, m, q, s in zip(pairs, ml, pid, sat):
            memo[frozenset(p)] = (float(m), float(q), bool(s))
        n_computed[0] += len(pairs)

    def bad(p):
        m, _, s = memo[frozenset(p)]
        return (not np.isfinite(m)) or s or m > a.dmax

    # ---------------------------------------------------------------- clade-closed selection
    t1 = time.time()
    state = {}
    for h in cand:
        ps = pres[h]
        cl = []
        for c in panel.internal:
            S = [s for s in ps if c in set(panel.anc[s])]
            if len(S) >= a.min_1to1:
                ag = panel.age[c]
                cl.append((-len(S), -(ag if np.isfinite(ag) else -1.0), c, S))
        cl.sort(key=lambda x: x[:3])
        state[h] = dict(cl=cl, ptr=0, k=a.k0, sel=None)
    if a.exhaustive:
        allp = []
        for h in cand:
            g1 = hog_stats[h]["one_gene"]
            allp += [(g1[x], g1[y]) for x, y in itertools.combinations(pres[h], 2)]
        compute(allp)
    pending = list(cand)
    rnd = 0
    while pending:
        rnd += 1
        req, nxt = [], []
        for h in pending:
            st = state[h]; g1 = hog_stats[h]["one_gene"]
            while st["ptr"] < len(st["cl"]):
                S = st["cl"][st["ptr"]][3]
                prs = [(x, y) for x, y in itertools.combinations(S, 2)]
                known = [p for p in prs if frozenset((g1[p[0]], g1[p[1]])) in memo]
                if any(bad((g1[x], g1[y])) for x, y in known):
                    st["ptr"] += 1; st["k"] = a.k0
                    continue
                unk = [p for p in prs if frozenset((g1[p[0]], g1[p[1]])) not in memo]
                if not unk:
                    st["sel"] = st["cl"][st["ptr"]]
                    break
                unk.sort(key=lambda p: -panel.pair_age(*p))
                req += [(g1[x], g1[y]) for x, y in unk[:st["k"]]]
                st["k"] *= 4
                nxt.append(h)
                break
        compute(req)
        print(f"selection round {rnd}: {len(pending)} pending, {len(req)} pairs requested, "
              f"{n_computed[0]} computed so far ({time.time()-t1:.0f}s)", flush=True)
        pending = nxt
    save_cache()
    rows = []
    for h in cand:
        st = state[h]; hs = hog_stats[h]; g1 = hs["one_gene"]
        r = dict(hog=h, n_genes=hs["n_genes"], n_genes_filtered=hs["n_genes_filtered"], n_species_1to1=len(hs["one"]),
                 presample=[taxid_of(s) for s in pres[h]], presample_genes=[g1[s] for s in pres[h]],
                 selected=[], selected_genes=[], sel_clade=None, root_age_mya=np.nan, stratum=pd.NA)
        if st["sel"] is not None:
            _, _, c, S = st["sel"]
            if len(S) > a.max_sel:
                S = greedy_cover(panel, S, a.max_sel, rng_for("hog_presample", crc(h), 1))
            m = panel.mrca(S)
            age = panel.age[m]
            r.update(selected=[taxid_of(s) for s in S], selected_genes=[g1[s] for s in S], sel_clade=panel.name[c],
                     sel_mrca=panel.name[m], root_age_mya=age,
                     stratum=next((k for k, lo, hi in STRATA if lo <= age < hi), pd.NA) if np.isfinite(age) else pd.NA,
                     mean_plddt=float(np.mean([gpl[g1[s]] for s in S])))
        rows.append(r)
    hogs = pd.DataFrame(rows)
    nsel = int(hogs.sel_clade.notna().sum())
    log.row("selection", step="clade_closed", hogs_in=len(cand), hogs_out=nsel,
            note=f"lazy={not a.exhaustive}; aa_ml_seq pairs computed {n_computed[0]}; rounds {rnd}")
    nst = int(hogs.stratum.notna().sum())
    log.row("selection", step="dated_root", hogs_in=nsel, hogs_out=nst)
    summary["presample_pairs_total"] = int(sum(len(v) * (len(v) - 1) // 2 for v in pres.values()))
    summary["seqonly_pairs_computed_selection"] = n_computed[0]

    # ---------------------------------------------------------------- labels + E6 helpers
    from ancestra.newick import TreeIndex
    from hog_pair_labels import label_pairs, parent_map
    tree = ox.species_tree
    tix, par = TreeIndex(tree), parent_map(tree)
    kids = {n.name: [c.name for c in n.children] for n in tree.preorder()}
    leaves_under = {n.name: set(n.leaf_names()) for n in tree.preorder()}
    is_leaf = {n.name: n.is_leaf() for n in tree.preorder()}
    topd = {g.id: g for g in tops}

    hidx = {}

    def lab_of(h, x, y):
        if h not in hidx:
            hidx.clear()           # keep one HOG index at a time (memory)
            hidx[h] = HogIndex(topd[h], sp_of, par)
        return hidx[h].label(x, y)

    def e6_candidates(h):
        """In-paralog pairs (non-leaf level) among filtered genes + O1/O2 sets; aa filters applied later."""
        gs = [x for x in hog_genes[h] if x in passed]
        bysp = defaultdict(list)
        for x in gs:
            bysp[sp_of[x]].append(x)
        single = {s: v[0] for s, v in bysp.items() if len(v) == 1}
        out = []
        for s, v in bysp.items():
            if len(v) < 2:
                continue
            for x, y in itertools.combinations(sorted(v), 2):
                rel, P, frm, _, dcg = lab_of(h, x, y)
                if rel != "in_paralog" or not P or is_leaf.get(P, True):
                    continue
                pP = par.get(P)
                if not pP:
                    continue
                sis1 = set().union(*[leaves_under[c] for c in kids[pP] if c != P]) if len(kids[pP]) > 1 else set()
                gp = par.get(pP)
                sis2 = set().union(*[leaves_under[c] for c in kids[gp] if c != pP]) if gp else set()

                def orth(sset):
                    o = []
                    for t in sorted(sset):
                        g = single.get(t)
                        if g and lab_of(h, x, g)[0] == "ortholog" and lab_of(h, y, g)[0] == "ortholog":
                            o.append(g)
                    return o
                O1, O2 = orth(sis1), orth(sis2)
                if O1:
                    out.append(dict(hog=h, species=s, a=x, b=y, P=P, frm=frm, dcg=dcg, O1=O1, O2=O2))
        return out

    # ---------------------------------------------------------------- sampling
    elig = hogs[hogs.stratum.notna()].copy()
    elig["stratum"] = elig.stratum.astype(int)
    perm = {}
    for k, _, _ in STRATA:
        hs = sorted(elig.loc[elig.stratum == k, "hog"])
        perm[k] = [hs[i] for i in rng_for("hog_sample", k).permutation(len(hs))]
    avail = {k: len(v) for k, v in perm.items()}

    def quotas(N):
        N = min(N, sum(avail.values()))
        base = N // 4
        q = {k: base for k in avail}
        extra = N - 4 * base
        for k in rng_for("hog_sample", 99).permutation(4)[:extra]:
            q[int(k) + 1] += 1
        while True:
            over = {k: q[k] - avail[k] for k in q if q[k] > avail[k]}
            if not over:
                return q
            left = sum(over.values())
            for k in over:
                q[k] = avail[k]
            spare = {k: avail[k] - q[k] for k in q if avail[k] > q[k]}
            tot = sum(spare.values())
            add = {k: left * v / tot for k, v in spare.items()}
            fl = {k: int(np.floor(v)) for k, v in add.items()}
            rem = left - sum(fl.values())
            fr = sorted(spare, key=lambda k: (-(add[k] - fl[k]), rng_for("hog_sample", 98, k).random()))
            for k in fr[:rem]:
                fl[k] += 1
            for k in fl:
                q[k] += fl[k]

    q = quotas(a.n_hogs)
    sample = [h for k in perm for h in perm[k][:q[k]]]
    for k in perm:
        log.row("sampling", step="initial_quota", stratum=k, eligible=avail[k], quota=q[k])
    # E6 + primary pair counts on the initial sample
    t1 = time.time()
    e6c = {h: e6_candidates(h) for h in sample}
    hidx.clear()
    print(f"labels + E6 candidates for {len(sample)} HOGs ({time.time()-t1:.0f}s)", flush=True)
    # aa filters for paralog pairs
    compute([(c["a"], c["b"]) for h in sample for c in e6c[h]])
    save_cache()
    e6n = Counter()
    e6_ok = {}
    for h in sample:
        keep = []
        for c in e6c[h]:
            m, pid, s = memo[frozenset((c["a"], c["b"]))]
            e6n["candidates"] += 1
            if not (np.isfinite(m) and not s and m <= a.dmax):
                e6n["fail_aa_ml_seq"] += 1; continue
            if pid > a.pid_max:
                e6n["fail_pid"] += 1; continue
            if a.e6_max_o:
                for key in ("O1", "O2"):
                    if len(c[key]) > a.e6_max_o:
                        e6n[f"capped_{key}"] += 1
                        r_ = rng_for("e6_p1_choice", crc(c["a"] + "|" + c["b"]), 1 if key == "O1" else 2)
                        c[key] = sorted(r_.choice(c[key], a.e6_max_o, replace=False).tolist())
            keep.append(c)
        if a.e6_max_per_dcg:
            byd = defaultdict(list)
            for c in keep:
                byd[c["dcg"]].append(c)
            keep = []
            for d_, cs in byd.items():
                if len(cs) > a.e6_max_per_dcg:
                    e6n["capped_dcg_triplets_dropped"] += len(cs) - a.e6_max_per_dcg
                    idx = rng_for("e6_p1_choice", crc(d_), 3).choice(len(cs), a.e6_max_per_dcg, replace=False)
                    cs = [cs[i] for i in sorted(idx)]
                keep += cs
        e6_ok[h] = keep
    hsel = hogs.set_index("hog")

    def npairs(h):
        n = len(hsel.at[h, "selected_genes"]); prim = n * (n - 1) // 2
        e6 = set()
        for c in e6_ok[h]:
            e6.add(frozenset((c["a"], c["b"])))
            for o in c["O1"] + c["O2"]:
                e6.add(frozenset((c["a"], o))); e6.add(frozenset((c["b"], o)))
        return prim, len(e6)
    pc = {h: npairs(h) for h in sample}
    tot = sum(p + e for p, e in pc.values())
    N = len(sample)
    summary["initial_sample"] = dict(n=N, primary=sum(p for p, _ in pc.values()), e6=sum(e for _, e in pc.values()))
    while tot > a.pair_budget and N > 0:
        N = int(np.floor(N * a.pair_budget / tot * 0.995))
        q = quotas(N)
        sample = [h for k in perm for h in perm[k][:q[k]]]
        tot = sum(sum(pc[h]) for h in sample)
        log.row("sampling", step="budget_rescale", n=N, total_pairs=tot)
    for k in perm:
        log.row("sampling", step="final_quota", stratum=k, eligible=avail[k], sampled=q[k])
    insample = set(sample)
    hogs["in_sample"] = hogs.hog.isin(insample)

    # ---------------------------------------------------------------- split (once)
    spl = {}
    for k in perm:
        hs = [h for h in perm[k][:q[k]]]
        p2 = rng_for("split", k).permutation(len(hs))
        nt = int(round(a.tune_frac * len(hs)))
        for r_, i in enumerate(p2):
            spl[hs[i]] = "tune" if r_ < nt else "confirm"
    if a.projection_only:
        print("projection only: split and pair lists not written")
    elif split_path.exists():
        old = pd.read_csv(split_path, sep="\t")
        od = dict(zip(old.hog, old.split))
        if od != spl:
            sys.exit(f"{split_path} exists and differs from the recomputed split ({len(od)} vs {len(spl)} HOGs) — "
                     "never regenerate; investigate")
        print(f"{split_path} exists and matches")
    else:
        st_of = dict(zip(hogs.hog, hogs.stratum))
        sdf = pd.DataFrame({"hog": list(spl), "stratum": [int(st_of[h]) for h in spl], "split": list(spl.values())})
        sdf.sort_values(["stratum", "hog"]).to_csv(split_path, sep="\t", index=False)
        import hashlib
        summary["split_md5"] = hashlib.md5(split_path.read_bytes()).hexdigest()
        print(f"wrote {split_path} md5 {summary['split_md5']}")
    hogs["split"] = [spl.get(h, "none") for h in hogs.hog]
    for k in perm:
        sub = [h for h in perm[k][:q[k]]]
        log.row("sampling", step="split", stratum=k, tune=sum(spl[h] == "tune" for h in sub),
                confirm=sum(spl[h] == "confirm" for h in sub))

    # ---------------------------------------------------------------- pair lists
    def meta(x, y):
        rel, lvl, frm, o2o, dcg = lab_of(cur_h[0], x, y)
        return rel, int(o2o), dcg, lvl, frm
    cur_h = [None]
    prim = []
    for h in sample:
        cur_h[0] = h
        gs = hsel.at[h, "selected_genes"]
        for x, y in itertools.combinations(sorted(gs), 2):
            rel, o2o, dcg, lvl, frm = meta(x, y)
            prim.append((x, y, "primary", h, spl[h], taxid_of(sp_of[x]), taxid_of(sp_of[y]), rel, o2o, dcg, lvl, frm))
    cols = ["gene_a", "gene_b", "pair_set", "hog", "split", "species_a", "species_b", "relation", "one_to_one",
            "dcg_id", "event_level", "dup_branch_from"]
    prim = pd.DataFrame(prim, columns=cols)
    e6rows, trip = {}, []
    for h in sample:
        cur_h[0] = h
        hidx.clear()
        for c in e6_ok[h]:
            x, y = (c["a"], c["b"]) if rng_for("e6_p1_choice", crc(c["a"] + "|" + c["b"])).random() < 0.5 \
                else (c["b"], c["a"])
            trip.append(dict(hog=h, split=spl[h], dcg_id=c["dcg"], species=taxid_of(c["species"]), P1=x, P2=y,
                             event_level=c["P"], dup_branch_from=c["frm"], O1=",".join(c["O1"]),
                             O2=",".join(c["O2"]), n_O1=len(c["O1"]), n_O2=len(c["O2"])))
            items = [((x, y), "e6_para")] + [((p, o), "e6_orth") for o in c["O1"] for p in (x, y)] + \
                    [((p, o), "e6_ldo") for o in c["O2"] for p in (x, y)]
            for (u, v), ps in items:
                key = frozenset((u, v))
                if key in e6rows:
                    continue
                rel, o2o, dcg, lvl, frm = meta(u, v)
                e6rows[key] = (u, v, ps, h, spl[h], taxid_of(sp_of[u]), taxid_of(sp_of[v]), rel, o2o, dcg, lvl, frm)
    e6 = pd.DataFrame(list(e6rows.values()), columns=cols)
    pk = set(map(frozenset, zip(prim.gene_a, prim.gene_b)))
    ov = e6.apply(lambda r: frozenset((r.gene_a, r.gene_b)) in pk, axis=1) if len(e6) else pd.Series([], dtype=bool)
    if ov.any():
        print(f"dropping {int(ov.sum())} E6 pairs already in the primary list")
        e6 = e6[~ov]
    trip = pd.DataFrame(trip)
    e6n["triplets"] = len(trip); e6n["triplets_with_O2"] = int((trip.n_O2 > 0).sum()) if len(trip) else 0
    e6n["dcg_groups"] = int(trip.dcg_id.nunique()) if len(trip) else 0
    e6n["hogs_with_triplets"] = int(trip.hog.nunique()) if len(trip) else 0
    for k, v in e6n.items():
        log.row("e6", step=k, n=v)
    for ps, n in e6.pair_set.value_counts().items():
        log.row("e6", step=f"pairs_{ps}", n=int(n))
    if not a.projection_only:
        pldir.mkdir(parents=True, exist_ok=True)
        prim.to_csv(pldir / "primary.tsv", sep="\t", index=False)
        e6.to_csv(pldir / "e6.tsv", sep="\t", index=False)
        trip.to_csv(pldir / "e6_triplets.tsv", sep="\t", index=False)
    # ---------------------------------------------------------------- hogs.parquet
    for c_ in ("misfit_ratio",):
        hogs[c_] = np.nan
    hogs["excl_misfit"] = pd.array([pd.NA] * len(hogs), dtype="boolean")
    hogs["cath_class"] = None          # filled by the §4.3 step
    for m in ("aa", "tdi", "lddt"):
        for al in ("tdiaa", "tmalign", "seqonly"):
            hogs[f"fail_frac_{m}_{al}"] = np.nan
            hogs[f"excl_{m}_{al}"] = pd.array([pd.NA] * len(hogs), dtype="boolean")
    hogs["stratum"] = hogs.stratum.astype("Int64")
    hogs["n_primary_pairs"] = [pc[h][0] if h in insample else 0 for h in hogs.hog]
    hogs["n_e6_pairs"] = [pc[h][1] if h in insample else 0 for h in hogs.hog]
    lead = ["hog", "n_species_1to1", "presample", "selected", "sel_clade", "root_age_mya", "stratum", "in_sample",
            "split", "cath_class", "mean_plddt", "misfit_ratio", "excl_misfit"]
    hogs = hogs[lead + [c_ for c_ in hogs.columns if c_ not in lead]]
    hogs.to_parquet(outdir / ("hogs_projection.parquet" if a.projection_only else "hogs.parquet"), index=False)
    by = hogs[hogs.in_sample].groupby("stratum")
    summary["final"] = dict(n_candidates=len(cand), n_selected=nsel, n_dated=nst, n_sampled=int(hogs.in_sample.sum()),
                            per_stratum_eligible=avail, per_stratum_sampled={int(k): int(v) for k, v in by.size().items()},
                            sel_size_mean=float(hogs.loc[hogs.in_sample, "selected"].map(len).mean()) if len(sample) else 0,
                            primary_pairs=len(prim), e6_pairs=len(e6), e6_triplets=len(trip), total_pairs=len(prim) + len(e6),
                            seqonly_pairs_computed_total=n_computed[0], wall_s=round(time.time() - t0, 1))
    summary["e6"] = dict(e6n)
    summary["eligible_root_age_by_stratum"] = {int(k): int(v) for k, v in elig.groupby("stratum").size().items()}
    json.dump(summary, open(outdir / ("hogs_projection_summary.json" if a.projection_only else "hogs_build_summary.json"),
                            "w"), indent=1, default=str)
    log.write()
    pool.close(); pool.join()
    print(json.dumps(summary["final"], indent=1, default=str))
    return 0


if __name__ == "__main__":
    sys.exit(main())
