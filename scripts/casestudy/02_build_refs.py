"""Curated Swiss-Prot reference set for subfamily labelling (X6 case study).
Sources (UniProt REST, 2026-10-30): reviewed entries of UniProt families
'actin family', 'FtsA/MreB family', 'tubulin family', 'FtsZ family', 'CetZ family',
plus ParM (P11904, reviewed) and BtubA/BtubB of Prosthecobacter dejongeii (Q8GCC5, Q8GCC1; unreviewed,
the only curated-literature BtubA/B entries). Labels are derived from UniProt subfamily
annotation first, then from the recommended protein name (rules below).
Output: refs/refs.fasta (>acc|label) and refs/refs_labels.tsv
"""
import pandas as pd, re, subprocess, io, requests, time
R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution/data/interim/casestudy/refs"

def lab_actin(name, fam):
    n = name.split(" (")[0]
    m = re.search(r"ARP(\d+) subfamily", fam)
    if "Plant ARP7" in fam: return "Arp7_plant"
    if m: return f"Arp{m.group(1)}"
    if "archaeal actin" in fam.lower() or "Archaeal actin" in n: return "Ta0583_archaeal_actin"
    if "Crenactin" in n: return "crenactin"
    if "POTE" in fam or "POTE" in n: return "actin_POTE"
    if re.search(r"subunit arp9", n): return "Arp9"
    if re.search(r"centractin", n, re.I): return "Arp1"
    m = re.search(r"Actin-related protein (\d+)", n)
    if m: return f"Arp{m.group(1)}"
    if re.search(r"Actin-like protein (6A|6B)|BAF53", n): return "Arp4"
    if re.search(r"Actin-like protein arp6|Actin-like protein ARP6", n, re.I): return "Arp6"
    if re.search(r"Actin-like protein (7A|7B|8|9|10)|Actin-related protein T|ACTRT|Actin-like protein 9", n): return "actin-like_metazoan_orphan"
    if re.search(r"Actin-related protein 3B", n): return "Arp3"
    if re.search(r"Actin-related protein 2", n): return "Arp2"
    if re.search(r"Actin-like", n): return "actin-like_other"
    return "actin"

def lab_mreb(name, gene):
    n = name.lower()
    if "ftsa" in n: return "FtsA"
    if "mamk" in n: return "MamK"
    return "MreB"

def lab_tub(name):
    n = name.lower()
    for k in ["alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta"]:
        if re.search(rf"\b{k}", n): return f"tubulin_{k}"
    return "tubulin_other"

rows = []
t = pd.read_csv(f"{R}/actin_fam.tsv", sep="\t").fillna("")
for _, r in t.iterrows(): rows.append((r.Entry, lab_actin(r["Protein names"], r["Protein families"]), r.Organism))
t = pd.read_csv(f"{R}/mreb_fam.tsv", sep="\t").fillna("")
for _, r in t.iterrows(): rows.append((r.Entry, lab_mreb(r["Protein names"], r["Gene Names (primary)"]), r.Organism))
t = pd.read_csv(f"{R}/tub_fam.tsv", sep="\t").fillna("")
for _, r in t.iterrows(): rows.append((r.Entry, lab_tub(r["Protein names"]), r.Organism))
t = pd.read_csv(f"{R}/ftsz_fam.tsv", sep="\t").fillna("")
for _, r in t.iterrows():
    f = r["Protein families"]
    rows.append((r.Entry, "TubZ" if "TubZ" in f else "PhuZ" if "PhuZ" in f else "FtsZ", r.Organism))
t = pd.read_csv(f"{R}/cetz_fam.tsv", sep="\t").fillna("")
for _, r in t.iterrows(): rows.append((r.Entry, "CetZ", r.Organism))
rows += [("P11904", "ParM", "Escherichia coli"), ("Q8GCC5", "BtubA", "Prosthecobacter dejongeii"),
         ("Q8GCC1", "BtubB", "Prosthecobacter dejongeii")]
lab = pd.DataFrame(rows, columns=["acc", "label", "organism"]).drop_duplicates("acc")
lab["superfamily"] = lab.label.map(lambda x: "tubulin" if x.startswith("tubulin") or x in
                                   ("FtsZ", "CetZ", "TubZ", "PhuZ", "BtubA", "BtubB") else "actin")
lab.to_csv(f"{R}/refs_labels.tsv", sep="\t", index=False)
print(lab.label.value_counts().to_string())
# sequences
seqs = {}
accs = list(lab.acc)
for i in range(0, len(accs), 100):
    for _ in range(5):
        r = requests.get("https://rest.uniprot.org/uniprotkb/accessions",
                         params={"accessions": ",".join(accs[i:i+100]), "format": "fasta", "size": 500}, timeout=120)
        if r.ok: break
        time.sleep(5)
    r.raise_for_status()
    acc = None
    for line in r.text.splitlines():
        if line.startswith(">"):
            acc = line.split("|")[1]; seqs[acc] = []
        else: seqs[acc].append(line)
L = dict(zip(lab.acc, lab.label))
with open(f"{R}/refs.fasta", "w") as fo:
    for a, s in seqs.items():
        if a in L: fo.write(f">{a}|{L[a]}\n{''.join(s)}\n")
print("refs with seq", len(seqs), "of", len(lab))
