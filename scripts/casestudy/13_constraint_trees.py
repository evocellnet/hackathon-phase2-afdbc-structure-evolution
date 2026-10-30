"""X6 node 2: constraint trees for the LBA sensitivity (constrained ML search, IQ-TREE -g).

actin  : ((eukaryotic actin+Arps),(Asgard actins)) as a clade, prokaryotic homologues free outside it
         -> the 'textbook' placement (Asgard actins sister to eukaryotic actin+Arps).
tubulin: eukaryotic tubulins (alpha..epsilon + divergent; organellar FtsZ excluded) form a clade that excludes
         the 3 artubulins and BtubA/B (free) -> artubulins/BtubAB cannot nest inside eukaryotic tubulins.
Writes data/interim/casestudy/constraint/{sf}_constraint.nwk
"""
import pandas as pd
R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
sel = pd.read_csv(f"{R}/results/casestudy/sample_members.tsv", sep="\t")
for sf in ("actin", "tubulin"):
    s = sel[sel.superfamily == sf]
    euk = s[(s.domain == "Eukaryota") & (s.subfamily != "euk_FtsZ")].seqid.tolist()
    if sf == "actin":
        asg = s[s.domain == "Asgard"].seqid.tolist()
        rest = s[~s.seqid.isin(euk + asg)].seqid.tolist()
        nwk = f"(({','.join(euk)}),({','.join(asg)})),{','.join(rest)}"
    else:
        rest = s[~s.seqid.isin(euk)].seqid.tolist()
        nwk = f"({','.join(euk)}),{','.join(rest)}"
    open(f"{R}/data/interim/casestudy/constraint/{sf}_constraint.nwk", "w").write(f"({nwk});\n")
    print(sf, len(euk), len(s))
