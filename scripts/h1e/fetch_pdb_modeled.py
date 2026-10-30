#!/usr/bin/env python
"""Per polymer chain: modeled (observed-in-coordinates) residue count and entity sequence length, via RCSB GraphQL.
SIFTS pdb_chain_uniprot segments span SEQRES incl. unobserved residues, so coverage must be corrected with this.
Entries: X-ray/EM with resolution <= 3.0 in pdb_entry_meta. -> data/raw/pdb_meta/pdb_chain_modeled.tsv.gz"""
import gzip, sys, time, requests, concurrent.futures as cf
R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution/data/raw/pdb_meta"
ids = []
with gzip.open(f"{R}/pdb_entry_meta.tsv.gz", "rt") as f:
    next(f)
    for l in f:
        p = l.rstrip("\n").split("\t")
        if p[2] in ("X-ray", "EM") and p[3] and float(p[3]) <= 3.0: ids.append(p[0].upper())
Q = ('query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id polymer_entities{entity_poly{rcsb_sample_sequence_length} '
     'polymer_entity_instances{rcsb_polymer_entity_instance_container_identifiers{auth_asym_id} rcsb_polymer_instance_info{modeled_residue_count}}}}}')
def run(chunk):
    for t in range(6):
        try:
            r = requests.post("https://data.rcsb.org/graphql", json={"query": Q, "variables": {"ids": chunk}}, timeout=180)
            r.raise_for_status(); d = r.json()
            if d.get("data") and d["data"].get("entries") is not None: return d["data"]["entries"]
        except Exception as e: err = e
        time.sleep(3 * (t + 1))
    raise RuntimeError(f"failed {chunk[:2]} {err}")
B = 100; n = 0
with gzip.open(f"{R}/pdb_chain_modeled.tsv.gz", "wt") as o, cf.ThreadPoolExecutor(8) as ex:
    o.write("pdb_id\tchain\tmodeled_residue_count\tentity_seq_length\n")
    for i, ents in enumerate(ex.map(run, [ids[i:i+B] for i in range(0, len(ids), B)])):
        for e in ents:
            if not e: continue
            for pe in e["polymer_entities"] or []:
                L = (pe["entity_poly"] or {}).get("rcsb_sample_sequence_length")
                for inst in pe["polymer_entity_instances"] or []:
                    ch = inst["rcsb_polymer_entity_instance_container_identifiers"]["auth_asym_id"]
                    mc = (inst["rcsb_polymer_instance_info"] or {}).get("modeled_residue_count")
                    o.write(f"{e['rcsb_id'].lower()}\t{ch}\t{'' if mc is None else mc}\t{'' if L is None else L}\n"); n += 1
        if i % 100 == 0: print(i, n, file=sys.stderr, flush=True)
print("entries", len(ids), "chains", n)
