"""X7 step 2: human kinome reference set for ePK group labels.
Source: UniProt REST (fetched 2026-10-30), query (organism_id:9606) AND (reviewed:true) AND
(family:"protein kinase superfamily"); fields accession, gene, SIMILARITY comment, DOMAIN features, sequence
-> data/interim/kinase/refs/human_pks_uniprot.tsv.
Group label = Manning-style group taken from the UniProt SIMILARITY family string:
  'AGC Ser/Thr' -> AGC, 'CAMK Ser/Thr' -> CAMK, 'CK1 Ser/Thr' -> CK1, 'CMGC Ser/Thr' -> CMGC, 'STE Ser/Thr' -> STE,
  'Tyr protein kinase' -> TK, 'TKL Ser/Thr' -> TKL; any other protein-kinase-superfamily family with a 'Protein
  kinase' DOMAIN feature (generic 'Ser/Thr protein kinase family', NEK, BUD32, PAN3) -> Other;
  'RIO-type' -> RIO, 'ADCK' -> ABC1 (atypical, kept to label eukaryotic atypicals); 'Alpha-type' -> excluded
  (no 'Protein kinase' domain feature / not PKL-Pkinase).
One reference sequence per 'Protein kinase' DOMAIN feature (kinase domain only; e.g. JAK pseudokinase + kinase domains
both kept, same label). RIO/ADCK proteins: UniProt has no 'Protein kinase' feature for some -> full sequence.
Output data/interim/kinase/refs/refs_kd.fasta (header acc_dN|label|gene) and refs_labels.tsv
"""
import re
import pandas as pd
R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
K = f"{R}/data/interim/kinase/refs"
d = pd.read_csv(f"{K}/human_pks_uniprot.tsv", sep="\t")
d.columns = ["acc", "gene", "sim", "dom", "seq"]
MAP = [("AGC Ser/Thr", "AGC"), ("CAMK Ser/Thr", "CAMK"), ("CK1 Ser/Thr", "CK1"), ("CMGC Ser/Thr", "CMGC"),
       ("STE Ser/Thr", "STE"), ("Tyr protein kinase", "TK"), ("TKL Ser/Thr", "TKL"), ("RIO-type", "RIO"),
       ("ADCK protein kinase", "ABC1"), ("Alpha-type", None)]


def label(s):
    s = str(s)
    m = re.search(r"elongs to the protein kinase superfamily\. ([^.]*)", s)
    fam = m.group(1) if m else ""
    for k, v in MAP:
        if k in fam:
            return v, fam
    return "Other", fam


rows, out = [], []
for r in d.itertuples():
    lab, fam = label(r.sim)
    if lab is None:
        continue
    doms = [(int(a), int(b)) for a, b in re.findall(r'DOMAIN (\d+)\.\.(\d+); /note="Protein kinase', str(r.dom))]
    if not doms:
        if lab in ("RIO", "ABC1"):
            doms = [(1, len(r.seq))]
        else:
            continue
    for i, (a, b) in enumerate(doms):
        out.append((f"{r.acc}_d{i+1}|{lab}|{r.gene}", r.seq[a - 1:b]))
        rows.append(dict(ref=f"{r.acc}_d{i+1}", acc=r.acc, gene=r.gene, label=lab, uniprot_family=fam, start=a, end=b))
with open(f"{K}/refs_kd.fasta", "w") as fo:
    for h, s in out:
        fo.write(f">{h}\n{s}\n")
t = pd.DataFrame(rows)
t.to_csv(f"{K}/refs_labels.tsv", sep="\t", index=False)
print(t.label.value_counts().to_string())
print("n refs", len(t), "proteins", t.acc.nunique())
