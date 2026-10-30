"""Kinase case study (X7), step 1: fetch all Pfam HMMs of clan CL0016 (protein kinase-like, PKL; 48 families,
Pfam 38.2 clan table in data/raw/pfam/Pfam-A.clans.tsv.gz) from the InterPro API and run
hmmsearch --cut_ga against all md5-matched panel proteins (180 euk + 100 prok; the 8 FASTA chunks built for X6 in
data/interim/casestudy/hmm_chunks/c*.fa = data/interim/casestudy/panel_all_v4match.fasta).
Outputs data/interim/kinase/pfam/{PFxxxxx.hmm, pkl.hmm}, data/interim/kinase/hmm/c*.domtbl
"""
import gzip, os, subprocess, urllib.request, hashlib
from concurrent.futures import ThreadPoolExecutor
R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
K = f"{R}/data/interim/kinase"
accs = [l.split("\t")[0] for l in gzip.open(f"{R}/data/raw/pfam/Pfam-A.clans.tsv.gz", "rt") if l.split("\t")[1] == "CL0016"]
print(len(accs), "CL0016 families")
import time
failed = []
for a in accs:
    p = f"{K}/pfam/{a}.hmm"
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        for attempt in range(4):
            try:
                raw = urllib.request.urlopen(f"https://www.ebi.ac.uk/interpro/api/entry/pfam/{a}?annotation=hmm",
                                             timeout=60).read()
                break
            except Exception as e:
                raw = None; time.sleep(5)
        if raw is None:
            failed.append(a); continue
        if raw[:2] == b"\x1f\x8b":
            raw = gzip.decompress(raw)
        open(p, "wb").write(raw)
print("HMM download failed (skipped, documented):", failed)
with open(f"{K}/pfam/pkl.hmm", "wb") as fo:
    for a in accs:
        if a in failed:
            continue
        fo.write(open(f"{K}/pfam/{a}.hmm", "rb").read())
print("pkl.hmm md5", hashlib.md5(open(f"{K}/pfam/pkl.hmm", "rb").read()).hexdigest())
subprocess.run(f"hmmpress -f {K}/pfam/pkl.hmm > /dev/null", shell=True, check=True)


def run(i):
    fa = f"{R}/data/interim/casestudy/hmm_chunks/c{i}.fa"
    out = f"{K}/hmm/c{i}"
    if os.path.exists(out + ".done"):
        return
    subprocess.run(f"hmmsearch --cut_ga --cpu 1 --noali --domtblout {out}.domtbl -o /dev/null {K}/pfam/pkl.hmm {fa}",
                   shell=True, check=True)
    open(out + ".done", "w").close()


with ThreadPoolExecutor(8) as ex:
    list(ex.map(run, range(8)))
print("done")
