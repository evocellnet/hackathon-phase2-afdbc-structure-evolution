#!/usr/bin/env python
"""Fetch single-chain PDB files (UniProt numbering) for H1e experimental controls.

usage: fetch_pdb_chains.py --acc P69905 P68871 ... | --acc-file list.txt
         [--table data/processed/h1e/pdb_chains.tsv] [--out DIR] [--workers 4]
For each UniProt accession: look up (pdb_id, auth chain) in the pdb_chains table, download the mmCIF (RCSB
files.rcsb.org/download/XXXX.cif.gz) and the per-residue SIFTS XML (EBI), keep chain residues that SIFTS maps to the
accession, renumber them to UniProt numbering (residue-level, via SIFTS PDBe/label_seq_id -> UniProt; robust to insertion
codes and engineered indels), keep first altloc, drop H/waters/ligands/unmapped residues (tags, linkers), rename MSE->MET,
write chain 'A'. Writes OUT/<acc>_<pdb>_<chain>.pdb (or .cif if a UniProt number > 9999) and OUT/fetch_summary.tsv with
sequence identity of the extracted residues against the UniProt FASTA (same positions)."""
import argparse, concurrent.futures as cf, gzip, os, sys, time, io, csv
import xml.etree.ElementTree as ET
import gemmi, requests
ROOT = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
NS = "{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}"

def http(url, tries=5):
    for k in range(tries):
        try:
            r = requests.get(url, timeout=120)
            if r.status_code == 404: return None
            r.raise_for_status(); return r.content
        except Exception:
            time.sleep(2 * (k + 1))
    raise RuntimeError("download failed: " + url)

def cached(url, path):
    if not os.path.exists(path):
        b = http(url)
        if b is None: raise FileNotFoundError(url)
        with open(path + ".tmp", "wb") as f: f.write(b)
        os.replace(path + ".tmp", path)
    return path

def sifts_map(xml_gz, chain, acc):
    """label_seq_id (PDBe resnum) -> UniProt resnum for residues of auth `chain` mapped to `acc`."""
    m = {}
    with gzip.open(xml_gz, "rb") as f:
        for _, el in ET.iterparse(f):
            if el.tag != NS + "residue": continue
            if el.get("dbSource") != "PDBe": el.clear(); continue
            ch = up = None
            for x in el.findall(NS + "crossRefDb"):
                s = x.get("dbSource")
                if s == "PDB": ch = x.get("dbChainId")
                elif s == "UniProt" and x.get("dbAccessionId") == acc: up = x.get("dbResNum")
            if ch == chain and up not in (None, "null"):
                try: m[int(el.get("dbResNum"))] = int(up)
                except ValueError: pass
            el.clear()
    return m

def uniprot_seq(acc, d):
    p = f"{d}/{acc}.fasta"
    if not os.path.exists(p):
        b = http(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta")
        open(p, "wb").write(b or b"")
    return "".join(l.strip() for l in open(p) if not l.startswith(">"))

def extract(acc, pdb, chain, out, cache):
    cif = cached(f"https://files.rcsb.org/download/{pdb.upper()}.cif.gz", f"{cache}/{pdb}.cif.gz")
    xml = cached(f"https://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/{pdb}.xml.gz", f"{cache}/{pdb}.sifts.xml.gz")
    mp = sifts_map(xml, chain, acc)
    st = gemmi.read_structure(cif); st.remove_hydrogens(); st.remove_alternative_conformations(); st.setup_entities()
    src = next((c for c in st[0] if c.name == chain), None)
    if src is None: raise ValueError("chain not in mmCIF")
    new = gemmi.Structure(); new.cell = st.cell; new.spacegroup_hm = st.spacegroup_hm; new.name = f"{acc}_{pdb}_{chain}"
    mod = gemmi.Model("1"); ch = gemmi.Chain("A"); seq = {}
    for r in src:
        if r.entity_type != gemmi.EntityType.Polymer or r.label_seq is None or r.label_seq not in mp: continue
        u = mp[r.label_seq]
        if u in seq: continue  # microheterogeneity duplicate
        nr = gemmi.Residue(); nr.name = "MET" if r.name == "MSE" else r.name; nr.seqid = gemmi.SeqId(u, " "); nr.het_flag = "A"
        for a in r: nr.add_atom(a)
        ch.add_residue(nr); seq[u] = gemmi.find_tabulated_residue(r.name).one_letter_code.upper() if gemmi.find_tabulated_residue(r.name) else "X"
    if not len(ch): raise ValueError("no mapped residues extracted")
    ch.residues.sort if False else None
    order = sorted(range(len(ch)), key=lambda i: ch[i].seqid.num)
    ch2 = gemmi.Chain("A")
    for i in order: ch2.add_residue(ch[i])
    mod.add_chain(ch2); new.add_model(mod); new.setup_entities()
    mx = max(seq)
    ext = "pdb" if mx <= 9999 else "cif"
    path = f"{out}/{acc}_{pdb}_{chain}.{ext}"
    if ext == "pdb": new.write_pdb(path)
    else: new.make_mmcif_document().write_file(path)
    return path, seq, mp

def one(a):
    acc, pdb, chain, out, cache = a
    try:
        path, seq, mp = extract(acc, pdb, chain, out, cache)
        us = uniprot_seq(acc, cache)
        match = sum(1 for u, aa in seq.items() if 1 <= u <= len(us) and us[u-1] == aa)
        nup = sorted(seq)
        return dict(uniprot=acc, pdb_id=pdb, chain=chain, n_extracted=len(seq), n_sifts_mapped=len(mp), n_match=match,
                    identity=round(match / len(seq), 4), uniprot_length=len(us), up_first=nup[0], up_last=nup[-1],
                    cov_extracted=round(len(seq) / len(us), 4) if us else "", file=path, status="ok")
    except Exception as e:
        return dict(uniprot=acc, pdb_id=pdb, chain=chain, status=f"FAIL: {e}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--acc", nargs="*", default=[]); ap.add_argument("--acc-file")
    ap.add_argument("--table", default=f"{ROOT}/data/processed/h1e/pdb_chains.tsv")
    ap.add_argument("--out", default=f"{ROOT}/data/interim/h1e/chains"); ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    accs = list(a.acc) + ([l.split()[0] for l in open(a.acc_file) if l.strip()] if a.acc_file else [])
    tab = {r["uniprot"]: r for r in csv.DictReader(open(a.table), delimiter="\t")}
    os.makedirs(a.out, exist_ok=True); cache = os.path.join(os.path.dirname(a.out), "cache"); os.makedirs(cache, exist_ok=True)
    jobs = []
    for x in accs:
        if x in tab: jobs.append((x, tab[x]["pdb_id"], tab[x]["chain"], a.out, cache))
        else: print(f"{x}: no qualifying chain in table", file=sys.stderr)
    with cf.ThreadPoolExecutor(a.workers) as ex: res = list(ex.map(one, jobs))
    keys = ["uniprot", "pdb_id", "chain", "n_extracted", "n_sifts_mapped", "n_match", "identity", "uniprot_length", "up_first", "up_last", "cov_extracted", "file", "status"]
    with open(f"{a.out}/fetch_summary.tsv", "w") as o:
        w = csv.DictWriter(o, keys, delimiter="\t", restval=""); w.writeheader(); w.writerows(res)
    for r in res: print("\t".join(str(r.get(k, "")) for k in keys))
