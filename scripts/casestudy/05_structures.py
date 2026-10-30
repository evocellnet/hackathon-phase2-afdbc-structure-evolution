"""X6 step 3a: extract AFDB v4 models of the sampled members from the panel FoldComp DBs, build a Foldseek DB
per superfamily, and write per-protein 3Di sequences (Foldseek 3Di, same residue indexing as the AA sequence,
md5 match guaranteed upstream)."""
import os, subprocess
import pandas as pd
import foldcomp

R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
CS = f"{R}/data/interim/casestudy"
sel = pd.read_csv(f"{R}/results/casestudy/sample_members.tsv", sep="\t")
dbs = {"Eukaryota": f"{R}/data/interim/structures/panel_v4"}
for sf in ["actin", "tubulin"]:
    d = sel[sel.superfamily == sf]
    od = f"{CS}/pdb/{sf}"
    os.makedirs(od, exist_ok=True)
    name2seqid = {f"AF-{a}-F1-model_v4": s for a, s in zip(d.acc, d.seqid)}
    got = set()
    for db, sub in [(f"{R}/data/interim/structures/panel_v4", d[d.domain == "Eukaryota"]),
                    (f"{R}/data/interim/structures/panel_prok_v4", d[d.domain != "Eukaryota"])]:
        ids = [f"AF-{a}-F1-model_v4" for a in sub.acc]
        with foldcomp.open(db, ids=ids) as fc:
            for name, pdb in fc:
                n = name.replace(".cif.gz", "").replace(".pdb", "")
                with open(f"{od}/{name2seqid[n]}.pdb", "w") as fo:
                    fo.write(pdb)
                got.add(n)
    miss = set(name2seqid) - got
    print(sf, "pdb written", len(got), "missing", len(miss), sorted(miss)[:5])
    fdb = f"{CS}/fsdb/{sf}"
    os.makedirs(f"{CS}/fsdb", exist_ok=True)
    subprocess.run(f"foldseek createdb {od} {fdb} --threads 4 -v 1", shell=True, check=True)
    subprocess.run(f"foldseek lndb {fdb}_h {fdb}_ss_h -v 1", shell=True, check=True)
    subprocess.run(f"foldseek convert2fasta {fdb}_ss {CS}/{sf}_3di.fasta -v 1", shell=True, check=True)
    subprocess.run(f"foldseek convert2fasta {fdb} {CS}/{sf}_aa_from_pdb.fasta -v 1", shell=True, check=True)
