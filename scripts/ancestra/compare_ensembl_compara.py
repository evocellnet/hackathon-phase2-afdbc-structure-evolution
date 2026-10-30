#!/usr/bin/env python3
"""External check of Ancestra pair labels against Ensembl Compara (gene-tree reconciliation).

For a random sample of human genes present in an Ancestra run (euk20), fetch Ensembl Compara
homologies (human paralogues + mouse orthologues/paralogues, with Compara's ``taxonomy_level``) via
the Ensembl REST API and compare:
  (1) human-mouse pairs: Ensembl ortholog_* vs Ancestra ortholog / out_paralog / not-in-same-HOG;
  (2) human-human paralog pairs: Ensembl duplication taxon vs Ancestra duplication level, both mapped
      to the guide-tree node directly below the duplication branch ("P", Ancestra convention): for an
      Ensembl duplication taxon T, P = the oldest guide-tree node on the human lineage inside T.

Usage: python scripts/ancestra/compare_ensembl_compara.py RUN_DIR N_GENES OUT_DIR
(RUN_DIR made by run_ancestra.py on UniProt reference proteomes incl. UP000005640_9606 and
UP000000589_10090.) Needs network (rest.ensembl.org, ~2 requests/gene).
"""
import gzip, json, os, random, re, sys, time, urllib.request
from collections import Counter
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ENS = os.path.join(ROOT, "data", "raw", "ancestra_test", "ensembl")
HUMAN, MOUSE = "UP000005640_9606", "UP000000589_10090"


def acc2ensg(path):
    m = {}
    with gzip.open(path, "rt") as fh:
        next(fh)
        for line in fh:
            acc, _, ens = line.rstrip("\n").split("\t")
            g = set(re.findall(r"(ENS(?:MUS)?G\d+)", ens))
            if len(g) == 1:
                m[acc] = g.pop()
    return m


def rest(url, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"Content-Type": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.load(r)
        except Exception as e:  # rate limit / transient
            time.sleep(2 + 3 * i)
    return None


def main():
    run, n, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    os.makedirs(out, exist_ok=True)
    h = acc2ensg(os.path.join(ENS, "UP000005640_ensembl.tsv.gz"))
    m = acc2ensg(os.path.join(ENS, "UP000000589_ensembl.tsv.gz"))
    g2acc_h, g2acc_m = {}, {}
    for a, g in h.items():
        g2acc_h.setdefault(g, []).append(a)
    for a, g in m.items():
        g2acc_m.setdefault(g, []).append(a)
    lab = pd.read_csv(os.path.join(run, "pair_labels.tsv.gz"), sep="\t", dtype=str)
    lab = lab[lab.species_a.isin([HUMAN, MOUSE]) & lab.species_b.isin([HUMAN, MOUSE])]
    key = {(r.gene_a, r.gene_b): r for r in lab.itertuples(index=False)}
    # sample human genes that are in the run and map to exactly one ENSG (one accession per ENSG)
    from ancestra.io import read_gene_table
    gt = read_gene_table(os.path.join(run, "ancestra", "stages", "input", "genes.tsv"))
    in_run = set(gt.names)
    cand = sorted(g for g, accs in g2acc_h.items() if len([a for a in accs if a in in_run]) == 1)
    random.Random(11).shuffle(cand)
    cand = cand[:n]
    cache = os.path.join(out, "ensembl_rest_cache.jsonl")
    done = {}
    if os.path.exists(cache):
        for line in open(cache):
            d = json.loads(line)
            done[d["gene"]] = d["data"]
    from concurrent.futures import ThreadPoolExecutor
    todo = [g for g in cand if g not in done]
    url = ("https://rest.ensembl.org/homology/id/human/{}?format=condensed;"
           "target_species=homo_sapiens;target_species=mus_musculus")
    with open(cache, "a") as fc, ThreadPoolExecutor(8) as ex:   # REST latency ~10-20 s/request
        for g, d in zip(todo, ex.map(lambda x: rest(url.format(x)), todo)):
            done[g] = d
            fc.write(json.dumps({"gene": g, "data": d}) + "\n")
            fc.flush()
    # taxonomy for level mapping
    from ancestra.newick import read_newick, TreeIndex
    from ancestra.speciestree import load_taxdump
    tax = load_taxdump(os.path.join(ROOT, "data", "interim", "ancestra", "taxdump"))
    parent, rank, name, merged = tax
    name2tid = {}
    for t, nm in name.items():
        name2tid.setdefault(nm, t)
    tree = read_newick(os.path.join(run, "guide_tree.nwk"))
    # human lineage in guide tree (leaf -> root) and NCBI taxid of each node
    node = {n.name: n for n in tree.preorder()}
    path = []
    x = node[HUMAN]
    while x is not None:
        path.append(x.name)
        x = x.parent
    def lineage(t):
        L = [t]
        while L[-1] != 1 and parent.get(L[-1], 1) != L[-1]:
            L.append(parent[L[-1]])
        return L
    node_tid = {HUMAN: 9606, MOUSE: 10090}
    for nm in path:
        if nm not in node_tid:
            base = re.sub(r"_\d+$", "", nm)
            node_tid[nm] = name2tid.get(base)
    hum_lin = set(lineage(9606))
    def ens_to_P(level_name):
        t = name2tid.get(level_name)
        if t is None:
            return None
        if level_name == "Homo sapiens":
            return HUMAN
        # oldest guide node on the human path whose taxon is inside T (T ancestor-or-equal of node)
        best = None
        for nm in path:          # leaf -> root
            tid = node_tid.get(nm)
            if tid is not None and t in lineage(tid):
                best = nm
        return best
    rows = []
    for g, d in done.items():
        if g not in cand or not d or not d.get("data"):
            continue
        accs_h = [a for a in g2acc_h.get(g, []) if a in in_run]
        if len(accs_h) != 1:
            continue
        a = accs_h[0]
        for hmo in d["data"][0]["homologies"]:
            tgt = hmo["id"]
            pool = g2acc_h if hmo["species"] == "homo_sapiens" else g2acc_m
            bs = [b for b in pool.get(tgt, []) if b in in_run]
            if len(bs) != 1:
                continue
            b = bs[0]
            k = (a, b) if a < b else (b, a)
            r = key.get(k)
            rows.append({"human_gene": g, "a": a, "b": b, "species": hmo["species"], "ens_type": hmo["type"],
                         "ens_level": hmo["taxonomy_level"], "ens_P": ens_to_P(hmo["taxonomy_level"]),
                         "anc_relation": r.relation if r is not None else "not_same_root_hog",
                         "anc_level": r.event_level if r is not None else ""})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(out, "ensembl_vs_ancestra_pairs.tsv"), sep="\t", index=False)
    rep = []
    rep.append(f"sampled human genes: {len(cand)}; homology pairs mapped to UniProt accessions in run: {len(df)}")
    hm = df[df.species == "mus_musculus"]
    rep.append("\n## human-mouse pairs: Ensembl type (rows) vs Ancestra relation (cols)\n")
    rep.append(pd.crosstab(hm.ens_type, hm.anc_relation, margins=True).to_string())
    hh = df[df.species == "homo_sapiens"]
    rep.append("\n## human-human Ensembl paralogues: Ancestra relation\n")
    rep.append(pd.crosstab(hh.ens_type, hh.anc_relation, margins=True).to_string())
    pp = hh[hh.anc_relation == "in_paralog"]
    agree = (pp.ens_P == pp.anc_level)
    rep.append(f"\n## duplication level (node below duplication branch) for pairs paralogous in both: "
               f"n={len(pp)}, exact agreement={agree.mean():.3f}")
    order = {nm: i for i, nm in enumerate(path)}   # 0 = human leaf (youngest)
    d_ = pp.assign(e=pp.ens_P.map(order), q=pp.anc_level.map(order)).dropna(subset=["e", "q"])
    rep.append(f"Ancestra younger than Ensembl: {(d_.q < d_.e).mean():.3f}; older: {(d_.q > d_.e).mean():.3f}; "
               f"within one node: {((d_.q - d_.e).abs() <= 1).mean():.3f}")
    rep.append("\nconfusion (rows Ensembl P, cols Ancestra P), human lineage order leaf->root:\n")
    ct = pd.crosstab(pd.Categorical(pp.ens_P, categories=path), pd.Categorical(pp.anc_level, categories=path))
    rep.append(ct.to_string())
    txt = "\n".join(rep)
    open(os.path.join(out, "ensembl_comparison.txt"), "w").write(txt + "\n")
    print(txt)


if __name__ == "__main__":
    main()
