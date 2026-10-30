#!/usr/bin/env python3
"""Panel species tree with TimeTree ages (ANALYSIS_PLAN §1.5, §2.2).

Topology
  1. Panel species whose NCBI taxid (strain or species) is a tip of the TimeTree-of-Life 2015 newick are mapped
     directly ("direct").
  2. Otherwise X = the deepest NCBI ancestor of the species that contains 2015 tips.
       - no other panel species in X                    -> any tip of X ("proxy_sole"; the induced placement is
                                                           then exact relative to the rest of the panel)
       - other panel species in X, a "neutral" tip exists (a tip of X in no sub-taxon of X that holds another
         panel species)                                 -> that tip ("proxy_group"); X is re-resolved in step 4
       - no neutral tip, or the species is undated against every other panel member of X
                                                        -> NCBI placement ("ncbi_attach"): polytomy child of the
                                                           MRCA of the panel members of the deepest NCBI taxon
                                                           holding another panel species
     Several panel species on one tip become a polytomy at that tip ("shared tip"), re-resolved in step 4.
  3. Induced subtree of the 2015 tree on the used tips.
  4. Re-resolution: for every proxy group X whose panel members form a clade in the tree and are dated pairwise,
     the clade is replaced by UPGMA on TimeTree 6 pairwise times (ties within 0.01 My -> polytomy). Inside such
     clades the 2015 placement of a proxy is arbitrary; TimeTree 6 is the dating source anyway.
Ages
  5. Each internal node: all spanning species pairs (one leaf below each of two different children).
     Representative pair = the dated spanning pair at the (lower) median of TimeTree 6 `median_age` (median of
     study estimates = the API median); lo95/hi95 = its CI. Nodes with no dated spanning pair -> NA.
  6. Monotonicity, top-down: a child older than (or as old as) its parent is set to the parent's age (logged).
  7. Collapse, top-down: an internal child whose internode (parent age - child age) < 5 My, or whose CI overlaps
     the parent's CI, is merged into the parent (logged). NA-aged nodes are not collapsed.
Outputs (TimeTree licence: everything with ages stays in data/interim/timetree/)
  panel_tree.nwk, node_ages.tsv, panel_tree_mapping.tsv, panel_tree_fixes.tsv, panel_tree_uncollapsed.nwk
  logs/filters/panel_tree.tsv (counts only)
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import treeutil as T  # noqa: E402
from taxo import Taxonomy  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
TTDIR = ROOT / "data/interim/timetree"
MIN_INTERNODE = 5.0
TIE_TOL = 0.01
FINE_RANKS = {"species", "subspecies", "strain", "species group", "species subgroup", "genus", "subgenus"}


def load_pairs(path):
    d = pd.read_csv(path, sep="\t")
    pa = {}
    for r in d.itertuples(index=False):
        lo, hi = float(r.ci_low), float(r.ci_high)
        if not (hi > lo):          # single-study pairs carry CI = (0, 0): no CI
            lo = hi = np.nan
        rec = None if pd.isna(r.age_mya) else (float(r.median_age), lo, hi,
                                               float(r.age_mya), float(r.adjusted_age), int(r.n_studies))
        pa[(r.leaf_a, r.leaf_b)] = rec
        pa[(r.leaf_b, r.leaf_a)] = rec
    return pa


def upgma(leaves, dist):
    """UPGMA with polytomies for ties (|d - dmin| <= TIE_TOL). Returns a Node."""
    cl = {i: (T.Node(l), [l]) for i, l in enumerate(leaves)}
    nxt = len(leaves)

    def d(a, b):
        return np.mean([dist(x, y) for x in cl[a][1] for y in cl[b][1]])
    D = {(a, b): d(a, b) for a, b in itertools.combinations(cl, 2)}
    while len(cl) > 1:
        dmin = min(D.values())
        # connected components of near-minimal links
        adj = defaultdict(set)
        for (a, b), v in D.items():
            if v <= dmin + TIE_TOL:
                adj[a].add(b); adj[b].add(a)
        start = min(adj)
        comp, st = set(), [start]
        while st:
            x = st.pop()
            if x not in comp:
                comp.add(x); st.extend(adj[x])
        n = T.Node()
        mem = []
        for c in sorted(comp):
            n.add(cl[c][0]); mem += cl[c][1]
            del cl[c]
        cl[nxt] = (n, mem)
        D = {k: v for k, v in D.items() if k[0] in cl and k[1] in cl}
        for c in cl:
            if c != nxt:
                D[(c, nxt)] = d(c, nxt)
        nxt += 1
    return next(iter(cl.values()))[0]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel", default=str(ROOT / "data/processed/panel/species_panel.tsv"))
    ap.add_argument("--tt2015", default=str(ROOT / "data/raw/timetree/TimetreeOfLife2015.nwk"))
    ap.add_argument("--tips", default=str(TTDIR / "ttol2015_tips_taxid.tsv"))
    ap.add_argument("--pairs", default=str(TTDIR / "pairwise_ages.tsv"))
    ap.add_argument("--outdir", default=str(TTDIR))
    ap.add_argument("--min-internode", type=float, default=MIN_INTERNODE)
    a = ap.parse_args(argv)
    out = Path(a.outdir)
    fixes = []                     # (step, node, detail)
    tx = Taxonomy()
    sp = pd.read_csv(a.panel, sep="\t")
    leaves = list(sp.leaf)
    lin = {r.leaf: tx.lineage(int(r.species_taxid)) for r in sp.itertuples()}
    pa = load_pairs(a.pairs)

    def dated(x, y):
        return pa.get((x, y)) is not None

    tt = pd.read_csv(a.tips, sep="\t").dropna(subset=["taxId"])
    tiplin, under = {}, defaultdict(list)
    for tip, tid in zip(tt.tip, tt.taxId.astype(int)):
        tiplin[tip] = tx.lineage(tid)
        for x in tiplin[tip]:
            under[x].append(tip)
    tip_of_taxid = {}
    for tip, l in tiplin.items():
        tip_of_taxid.setdefault(l[0], tip)
    # ---- 1/2 mapping
    mp = {}
    for r in sp.itertuples():
        for t in (int(r.taxid), int(r.species_taxid)):
            t = tx.norm(t)
            if t in tip_of_taxid:
                mp[r.leaf] = dict(mode="direct", tip=tip_of_taxid[t], taxon=t)
                break
    used = {m["tip"] for m in mp.values()}
    panel_in = defaultdict(set)    # NCBI taxon -> panel leaves
    for l, L in lin.items():
        for x in L:
            panel_in[x].add(l)
    for r in sp.itertuples():
        l = r.leaf
        if l in mp:
            continue
        L = lin[l]
        X = next(x for x in L if x in under)
        others = panel_in[X] - {l}
        if not others:
            tip = sorted(t for t in under[X] if t not in used)[0] if any(t not in used for t in under[X]) \
                else sorted(under[X])[0]
            mp[l] = dict(mode="proxy_sole", tip=tip, taxon=X)
        else:
            undated_all = all(not dated(l, o) for o in others)
            below = L[:L.index(X)]                      # p's own lineage below X
            occupied = set()
            for o in others:
                Lo = lin[o]
                occupied |= set(Lo[:Lo.index(X)]) if X in Lo else set()
            neutral = sorted(t for t in under[X] if not (set(tiplin[t][:tiplin[t].index(X)]) & occupied)
                             and t not in used)
            # prefer neutral tips sharing p's own sub-taxa
            neutral.sort(key=lambda t: -len(set(tiplin[t]) & set(below)))
            shared = [o for o in others if o in mp and mp[o]["mode"] == "direct"
                      and tx.rank.get(X, "") in FINE_RANKS]
            if undated_all:
                mp[l] = dict(mode="ncbi_attach", tip=None, taxon=X)
            elif neutral:
                mp[l] = dict(mode="proxy_group", tip=neutral[0], taxon=X)
            elif tx.rank.get(X, "") in FINE_RANKS or len(under[X]) <= 3:
                # congeners with too few tips: share the tip of a panel member (re-resolved below)
                tgt = sorted(shared)[0] if shared else None
                tip = mp[tgt]["tip"] if tgt else sorted(under[X])[0]
                mp[l] = dict(mode="shared_tip", tip=tip, taxon=X)
            else:
                mp[l] = dict(mode="ncbi_attach", tip=None, taxon=X)
        if mp[l]["tip"]:
            used.add(mp[l]["tip"])
    # ---- 3 induced tree
    tt15 = T.parse_newick(Path(a.tt2015).read_text())
    usedtips = {m["tip"] for m in mp.values() if m["tip"]}
    tr = T.prune_to(tt15, usedtips)
    bytip = defaultdict(list)
    for l, m in mp.items():
        if m["tip"]:
            bytip[m["tip"]].append(l)
    for n in list(tr.preorder()):
        if n.is_leaf():
            ls = sorted(bytip[n.name])
            if len(ls) == 1:
                n.name = ls[0]
            else:
                n.name = ""
                for x in ls:
                    n.add(T.Node(x))
    for n in tr.preorder():
        n.length = None
    # ncbi attachments
    for l, m in sorted(mp.items()):
        if m["mode"] != "ncbi_attach":
            continue
        L = lin[l]
        placed = {x.name for x in tr.leaves()}
        X2 = next(x for x in L[1:] if panel_in[x] & placed - {l})
        nodes, ls = T.index_tree(tr)
        M = T.mrca(nodes, ls, panel_in[X2] & placed)
        if M.is_leaf():          # single panel relative: make a cherry
            par = M.parent
            nn = T.Node(); i = par.children.index(M); par.children[i] = nn; nn.parent = par
            nn.add(M)
            M = nn
        M.add(T.Node(l))
        m["attach_taxon"] = X2
        fixes.append(("ncbi_attach", l, f"attached as polytomy child of MRCA({tx.name.get(X2, X2)}; "
                                        f"{len(panel_in[X2] & placed)} panel species)"))
    _lv = sorted(x.name for x in tr.leaves())
    assert _lv == sorted(leaves), ("leaf set mismatch", set(_lv) ^ set(leaves), len(_lv))
    # ---- 4 re-resolution of proxy groups
    groups = sorted({m["taxon"] for m in mp.values() if m["mode"] in ("proxy_group", "shared_tip")}
                    | {m["attach_taxon"] for m in mp.values() if m["mode"] == "ncbi_attach"},
                    key=lambda X: len(panel_in[X]))
    for X in groups:
        mem = panel_in[X]
        nodes, ls = T.index_tree(tr)
        M = T.mrca(nodes, ls, mem)
        allpairs = list(itertools.combinations(sorted(mem), 2))
        if ls[M.id] != frozenset(mem):
            fixes.append(("reresolve_skip", tx.name.get(X, X), f"not monophyletic in tree ({len(ls[M.id])} vs {len(mem)})"))
            continue
        if not all(dated(x, y) for x, y in allpairs):
            fixes.append(("reresolve_skip", tx.name.get(X, X), "undated pairs inside; 2015/proxy structure kept"))
            continue
        new = upgma(sorted(mem), lambda x, y: pa[(x, y)][3])
        par = M.parent
        if par is None:
            tr = new
        else:
            i = par.children.index(M)
            par.children[i] = new; new.parent = par
        fixes.append(("reresolve", tx.name.get(X, X), f"{len(mem)} species re-resolved by UPGMA on TimeTree 6 ages"))
    tr.parent = None
    T.prune_to  # noqa
    # suppress unary nodes (none expected)
    uncollapsed_n = sum(1 for n in tr.preorder() if not n.is_leaf())

    # ---- 5 ages
    def node_age(n, ls):
        kids = [ls[c.id] for c in n.children]
        recs, und = [], 0
        for i, j in itertools.combinations(range(len(kids)), 2):
            for x in kids[i]:
                for y in kids[j]:
                    r = pa.get((x, y))
                    if r is None:
                        und += 1
                    else:
                        recs.append((r[0], x, y, r))
        if not recs:
            return None, und, 0
        recs.sort()
        rep = recs[(len(recs) - 1) // 2]
        return dict(median=rep[0], lo95=rep[3][1], hi95=rep[3][2], adjusted=rep[3][4], n_studies=rep[3][5],
                    rep_pair=f"{rep[1]}|{rep[2]}", span_min=recs[0][0], span_max=recs[-1][0]), und, len(recs)

    def date_all(tr):
        nodes, ls = T.index_tree(tr)
        info = {}
        for n in nodes:
            if n.is_leaf():
                n.age = 0.0
                info[n.id] = dict(median=0.0, lo95=0.0, hi95=0.0, source="tip")
                continue
            r, und, nd = node_age(n, ls)
            if r is None:
                n.age = None
                info[n.id] = dict(median=np.nan, lo95=np.nan, hi95=np.nan, source="NA:no_dated_spanning_pair",
                                  n_span_dated=0, n_span_undated=und)
            else:
                n.age = r["median"]
                info[n.id] = dict(r, source=f"timetree6_api_median:{r['rep_pair']}", n_span_dated=nd,
                                  n_span_undated=und)
        return nodes, ls, info

    def ncbi_lca_name(lset):
        Ls = [lin[x] for x in lset]
        common = set(Ls[0]).intersection(*map(set, Ls[1:]))
        best = next(x for x in Ls[0] if x in common)
        return tx.name.get(best, str(best))

    nodes, ls, info = date_all(tr)
    # ---- 6 monotone top-down
    for n in tr.preorder():
        if n.is_leaf() or n.parent is None or n.age is None:
            continue
        p = n.parent
        while p is not None and p.age is None:
            p = p.parent
        if p is not None and n.age >= p.age:
            fixes.append(("monotone", ncbi_lca_name(ls[n.id]),
                          f"age {n.age:.2f} >= ancestor {ncbi_lca_name(ls[p.id])} {p.age:.2f}; set to {p.age:.2f}"))
            info[n.id]["source"] += f";monotone_fix_from_{n.age:.2f}"
            n.age = p.age
            info[n.id]["median"] = p.age
            info[n.id]["lo95"] = np.fmin(info[n.id]["lo95"], p.age); info[n.id]["hi95"] = np.fmax(info[n.id]["hi95"], p.age)
    (out / "panel_tree_uncollapsed.nwk").write_text(T.to_newick(tr, lengths=False, internal_names=False) + "\n")
    # ---- 7 collapse top-down
    n_coll = 0
    st = [tr]
    while st:
        p = st.pop()
        changed = True
        while changed:
            changed = False
            for c in list(p.children):
                if c.is_leaf() or c.age is None or p.age is None:
                    continue
                ci_c, ci_p = info[c.id], info[p.id]
                short = (p.age - c.age) < a.min_internode
                ovl = (ci_c["hi95"] >= ci_p["lo95"]) and not (np.isnan(ci_c["hi95"]) or np.isnan(ci_p["lo95"]))
                if short or ovl:
                    why = "internode<%.0fMy" % a.min_internode if short else "CI_overlap"
                    fixes.append(("collapse", ncbi_lca_name(ls[c.id]),
                                  f"{why}: child {c.age:.1f} [{ci_c['lo95']:.1f},{ci_c['hi95']:.1f}] into parent "
                                  f"{ncbi_lca_name(ls[p.id])} {p.age:.1f} [{ci_p['lo95']:.1f},{ci_p['hi95']:.1f}]"))
                    i = p.children.index(c)
                    for g in c.children:
                        g.parent = p
                    p.children[i:i + 1] = c.children
                    n_coll += 1
                    changed = True
        st.extend(c for c in p.children if not c.is_leaf())
    # final ids / names / lengths
    old = {id(n): info[n.id] for n in tr.preorder()}
    nodes, ls = T.index_tree(tr)
    rows = []
    k = 0
    for n in nodes:
        r = old[id(n)]
        if not n.is_leaf():
            k += 1
            n.name = f"N{k:03d}"
        nm = n.name
        rows.append(dict(node=nm, median=r["median"], lo95=r["lo95"], hi95=r["hi95"], source=r["source"],
                         is_leaf=n.is_leaf(), parent=n.parent.name if n.parent else "",
                         n_leaves=len(ls[n.id]), n_children=len(n.children),
                         ncbi_lca=ncbi_lca_name(ls[n.id]), adjusted_age=r.get("adjusted", np.nan),
                         n_studies=r.get("n_studies", np.nan), n_span_dated=r.get("n_span_dated", np.nan),
                         n_span_undated=r.get("n_span_undated", np.nan), span_min=r.get("span_min", np.nan),
                         span_max=r.get("span_max", np.nan),
                         leaves=",".join(sorted(ls[n.id])) if not n.is_leaf() else n.name))
    for n in nodes:
        if n.parent is not None and n.age is not None and n.parent.age is not None:
            n.length = n.parent.age - n.age
        else:
            n.length = None
    na = pd.DataFrame(rows)
    # leaves with an NCBI-attachment (no own dated split) are flagged
    att = {l for l, m in mp.items() if m["mode"] == "ncbi_attach"}
    na["attach"] = [("ncbi_polytomy" if x in att else "") for x in na.node]
    na.to_csv(out / "node_ages.tsv", sep="\t", index=False, float_format="%.4f")
    (out / "panel_tree.nwk").write_text(T.to_newick(tr) + "\n")
    mt = []
    for r in sp.itertuples():
        m = mp[r.leaf]
        mt.append(dict(leaf=r.leaf, taxid=r.taxid, species_taxid=r.species_taxid, name=r.name, group=r.group,
                       mode=m["mode"], tt2015_tip=m["tip"] or "", proxy_taxon=tx.name.get(m["taxon"], ""),
                       proxy_rank=tx.rank.get(m["taxon"], ""),
                       attach_taxon=tx.name.get(m.get("attach_taxon"), "") if m.get("attach_taxon") else ""))
    pd.DataFrame(mt).to_csv(out / "panel_tree_mapping.tsv", sep="\t", index=False)
    pd.DataFrame(fixes, columns=["step", "node", "detail"]).to_csv(out / "panel_tree_fixes.tsv", sep="\t", index=False)
    # ---- checks
    nodes, ls = T.index_tree(tr)
    ok_mono = all(n.parent is None or n.age is None or n.parent.age is None or n.age < n.parent.age
                  for n in nodes if not n.is_leaf())
    hm = T.mrca(nodes, ls, {"UP000005640_9606", "UP000000589_10090"})
    hmr = na.set_index("node").loc[hm.name]
    summ = dict(n_leaves=len(ls[0]), n_internal_uncollapsed=uncollapsed_n,
                n_internal=sum(1 for n in nodes if not n.is_leaf()), n_collapsed=n_coll,
                n_na_nodes=int((na.source.str.startswith("NA")).sum()),
                modes=pd.DataFrame(mt)["mode"].value_counts().to_dict(),
                n_monotone_fixes=sum(1 for f in fixes if f[0] == "monotone"),
                n_reresolved=sum(1 for f in fixes if f[0] == "reresolve"), monotone_ok=bool(ok_mono),
                human_mouse=dict(node=hm.name, median=float(hmr["median"]), lo95=float(hmr.lo95), hi95=float(hmr.hi95)),
                human_mouse_check_87pm10=bool(abs(hmr["median"] - 87) <= 10))
    (ROOT / "logs/filters").mkdir(parents=True, exist_ok=True)
    pd.DataFrame([(k, json.dumps(v) if isinstance(v, dict) else v) for k, v in summ.items() if k != "human_mouse"],
                 columns=["item", "value"]).to_csv(ROOT / "logs/filters/panel_tree.tsv", sep="\t", index=False)
    json.dump(summ, open(out / "panel_tree_summary.json", "w"), indent=1)
    print(json.dumps(summ, indent=1))
    assert ok_mono
    return 0


if __name__ == "__main__":
    sys.exit(main())
