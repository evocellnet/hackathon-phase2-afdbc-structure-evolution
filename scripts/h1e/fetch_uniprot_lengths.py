#!/usr/bin/env python
"""UniProt accession/length/organism for all entries with structure_3d:true in the H1e taxa
(paged UniProt REST search; taxonomy_id includes descendant strains). -> data/raw/pdb_meta/uniprot_3d_<taxid>.tsv"""
import sys, re, time, requests, concurrent.futures as cf
OUT = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution/data/raw/pdb_meta"
# usage: [all|genes]; "genes" adds gene_primary and writes uniprot_3d_genes_<taxid>.tsv for the vertebrate taxa (symbol-based pair forecast)
GENES = len(sys.argv) > 1 and sys.argv[1] == "genes"
FIELDS = "accession,length,organism_id,reviewed" + (",gene_primary" if GENES else "")
TAX = [9606, 10090, 10116, 559292, 3702, 7227, 6239, 7955, 83333, 9031, 9913]
def get(url, params):
    for k in range(6):
        try:
            r = requests.get(url, params=params, timeout=120); r.raise_for_status(); return r
        except Exception: time.sleep(3 * (k + 1))
    raise RuntimeError(url)
def run(t):
    q = "(structure_3d:true)" if t == "all" else f"(taxonomy_id:{t}) AND (structure_3d:true)"
    url = "https://rest.uniprot.org/uniprotkb/search"
    p = {"query": q, "fields": FIELDS, "format": "tsv", "size": 500}
    n = 0; first = True
    with open(f"{OUT}/uniprot_3d_{'genes_' if GENES else ''}{t}.tsv", "w") as o:
        while url:
            r = get(url, p if first else None)
            lines = r.text.splitlines()
            o.write("\n".join(lines if first else lines[1:]) + "\n"); n += len(lines) - 1; first = False
            m = re.search(r'<([^>]+)>; rel="next"', r.headers.get("Link", "")); url = m.group(1) if m else None
            print(t, n, file=sys.stderr, flush=True)
    return t, n
# usage: python fetch_uniprot_lengths.py [all]   ("all" = every UniProt entry with a 3D structure, any taxon)
todo = ["all"] if sys.argv[1:] == ["all"] else ([9606, 10090, 10116, 9913, 9031, 7955] if GENES else TAX)
with cf.ThreadPoolExecutor(6) as ex:
    for t, n in ex.map(run, todo): print(t, n)
