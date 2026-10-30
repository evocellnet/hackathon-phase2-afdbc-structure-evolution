"""X7 step 5: extract AFDB v4 models of the sampled proteins from the panel FoldComp DBs, crop each to its kinase
domain envelope (env_from..env_to, original residue numbering kept), build a Foldseek DB of the cropped domains and
write per-domain AA and 3Di sequences (3Di computed on the cropped domain, so the same structure input is used for
the projected-3Di alignments and for FoldMason).
Outputs data/interim/kinase/pdb/<seqid>.pdb, fsdb/pkl, pkl_3di.fasta, pkl_aa_from_pdb.fasta
"""
import os, subprocess
import pandas as pd
import foldcomp

R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
K = f"{R}/data/interim/kinase"
s = pd.read_csv(f"{R}/results/kinase/sample_members.tsv", sep="\t")
od = f"{K}/pdb"
os.makedirs(od, exist_ok=True)
info = {f"AF-{r.acc}-F1-model_v4": r for r in s.itertuples()}
got = set()
for db, sub in [(f"{R}/data/interim/structures/panel_v4", s[s.domain == "Eukaryota"]),
                (f"{R}/data/interim/structures/panel_prok_v4", s[s.domain != "Eukaryota"])]:
    ids = [f"AF-{a}-F1-model_v4" for a in sub.acc]
    with foldcomp.open(db, ids=ids) as fc:
        for name, pdb in fc:
            n = name.replace(".cif.gz", "").replace(".pdb", "")
            r = info[n]
            out = []
            for line in pdb.splitlines():
                if line.startswith(("ATOM", "HETATM")):
                    resi = int(line[22:26])
                    if r.env_from <= resi <= r.env_to:
                        out.append(line)
            out.append("END")
            with open(f"{od}/{r.seqid}.pdb", "w") as fo:
                fo.write("\n".join(out) + "\n")
            got.add(n)
miss = set(info) - got
print("pdb written", len(got), "missing", len(miss), sorted(miss)[:5])
assert not miss
os.makedirs(f"{K}/fsdb", exist_ok=True)
fdb = f"{K}/fsdb/pkl"
subprocess.run(f"foldseek createdb {od} {fdb} --threads 8 -v 1", shell=True, check=True)
subprocess.run(f"foldseek lndb {fdb}_h {fdb}_ss_h -v 1", shell=True, check=True)
subprocess.run(f"foldseek convert2fasta {fdb}_ss {K}/pkl_3di.fasta -v 1", shell=True, check=True)
subprocess.run(f"foldseek convert2fasta {fdb} {K}/pkl_aa_from_pdb.fasta -v 1", shell=True, check=True)
