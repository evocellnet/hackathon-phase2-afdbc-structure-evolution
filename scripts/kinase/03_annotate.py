"""X7 step 3: PKL membership, kinase-domain boundaries, domain pLDDT and ePK group labels.

Inputs: data/interim/kinase/hmm/c*.domtbl (hmmsearch --cut_ga, 48 CL0016 HMMs, all md5-matched panel proteins),
        data/interim/structures/panel{,_prok}_v4_plddt.npz (per-residue pLDDT of the AFDB v4 models),
        data/interim/kinase/refs/refs_kd.fasta (human kinase domains, 02_refs.py).
Rules
  * PKL member: >= 1 GA domain hit to any CL0016 family. The protein's kinase domain = the single highest-scoring
    CL0016 domain hit (clan competition: overlapping hits of other CL0016 families are ignored); its family is the
    protein's 'pfam_family'. Boundaries = envelope (env_from..env_to). n_pkl_domains = number of non-overlapping
    CL0016 domain hits (multi-kinase proteins, e.g. JAK, RSK, are kept but only the best domain is used).
  * 'ePK' = best family Pkinase (PF00069) or PK_Tyr_Ser-Thr (PF07714). Atypical/ELK families sampled for the tree:
    RIO1, ABC1, APH, Choline_kinase, Kdo (others are recorded but not sampled).
  * quality ('pass'): HMM coverage of the best hit >= 0.70 of model length; domain mean pLDDT >= 70 (mean of
    per-residue pLDDT over env_from..env_to); domain length <= 2x the family median of HMM-covered hits.
  * ePK group: DIAMOND blastp --ultra-sensitive best hit (bitscore) of the DOMAIN sequence vs the human kinase-domain
    references; label AGC/CAMK/CK1/CMGC/STE/TK/TKL/Other (UniProt family string; 02_refs.py). Second-best distinct
    label and bitscore margin recorded. Eukaryotic best-hit-TK domains outside Holozoa (Metazoa + Holozoa_other) are
    labelled 'ePK_TKlike_nonholozoan' (median identity 30 % vs 74 % in Metazoa; true PTKs are holozoan-specific).
    Note: eukaryotic 'Kdo_euk' members are all TP53RK/Bud32-like (best human hit TP53RK), i.e. the Bud32 atypical
    kinase, which the Kdo HMM out-scores. Prokaryotic ePKs keep 'prok_ePK' as class (best human label recorded only).
Outputs results/kinase/membership_all.tsv.gz, results/kinase/membership_summary.tsv,
        data/interim/kinase/pkl_domains.fasta (domain sequences of 'pass' proteins)
"""
import glob, subprocess
import numpy as np, pandas as pd

R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
K = f"{R}/data/interim/kinase"
OUT = f"{R}/results/kinase"
EPK = {"Pkinase", "PK_Tyr_Ser-Thr"}
ELK = {"RIO1", "ABC1", "APH", "Choline_kinase", "Kdo"}

rows = []
for f in sorted(glob.glob(f"{K}/hmm/c*.domtbl")):
    for line in open(f):
        if line.startswith("#"):
            continue
        x = line.split()
        rows.append((x[0], x[3], x[4].split(".")[0], int(x[2]), int(x[5]), float(x[13]), int(x[15]), int(x[16]),
                     int(x[19]), int(x[20])))
dom = pd.DataFrame(rows, columns=["target", "pfam_family", "pfam_acc", "tlen", "qlen", "dom_score", "hmm_from",
                                  "hmm_to", "env_from", "env_to"])
dom.to_csv(f"{K}/hmm_domtbl_all.tsv.gz", sep="\t", index=False)
print("domain hits", len(dom), "proteins", dom.target.nunique())


def n_nonoverlap(g):
    g = g.sort_values("dom_score", ascending=False)
    kept = []
    for a, b in zip(g.env_from, g.env_to):
        if all(b < c or a > d for c, d in kept):
            kept.append((a, b))
    return len(kept)


nd = dom.groupby("target").apply(n_nonoverlap)
best = dom.sort_values("dom_score", ascending=False).groupby("target").head(1).set_index("target")
m = best.reset_index()
m["n_pkl_domains"] = m.target.map(nd)
m["arch_pkl"] = m.target.map(dom.sort_values("env_from").groupby("target").pfam_family.apply(
    lambda s: "+".join(dict.fromkeys(s))))
m["hmm_cov"] = (m.hmm_to - m.hmm_from + 1) / m.qlen
m["dom_len"] = m.env_to - m.env_from + 1
m["acc"] = m.target.str.split("|").str[0]
m["leaf"] = m.target.str.split("|").str[1]

pt = pd.concat([pd.read_parquet(f"{R}/data/processed/panel/{f}", columns=["acc", "leaf", "len_v4", "mean_plddt"])
                for f in ["protein_table.parquet", "protein_table_prok.parquet"]])
sp = pd.concat([pd.read_csv(f"{R}/data/processed/panel/{f}", sep="\t")
                for f in ["species_panel.tsv", "species_panel_prok.tsv"]])[["leaf", "taxid", "name", "superkingdom",
                                                                               "group"]]
m = m.merge(pt, on=["acc", "leaf"], how="left").merge(sp, on="leaf", how="left")
m["domain"] = m.superkingdom
m.loc[m.group == "Promethearchaeota", "domain"] = "Asgard"

# per-residue pLDDT over the domain envelope
m["dom_plddt"] = np.nan
for tag, sk in (("panel_v4", "Eukaryota"), ("panel_prok_v4", None)):
    z = np.load(f"{R}/data/interim/structures/{tag}_plddt.npz")
    idx = {a: i for i, a in enumerate(z["acc"])}
    off, pl = z["offsets"], z["plddt"]
    sel = (m.superkingdom == "Eukaryota") if sk else (m.superkingdom != "Eukaryota")
    vals = []
    for r in m[sel].itertuples():
        i = idx.get(r.acc)
        if i is None:
            vals.append(np.nan); continue
        v = pl[off[i]:off[i + 1]].astype(float)
        assert len(v) == r.tlen, (r.acc, len(v), r.tlen)
        vals.append(v[r.env_from - 1:r.env_to].mean())
    m.loc[sel, "dom_plddt"] = vals
med = m[m.hmm_cov >= 0.7].groupby("pfam_family").dom_len.median()
m["cov_ok"] = m.hmm_cov >= 0.70
m["plddt_ok"] = m.dom_plddt >= 70
m["len_ok"] = m.dom_len <= 2 * m.pfam_family.map(med)
m["pass"] = m.cov_ok & m.plddt_ok & m.len_ok

# domain sequences
want = set(m.target)
seqs, cur = {}, None
for line in open(f"{R}/data/interim/casestudy/panel_all_v4match.fasta"):
    if line.startswith(">"):
        cur = line[1:].strip()
        cur = cur if cur in want else None
    elif cur:
        seqs[cur] = seqs.get(cur, "") + line.strip()
cand = m[m.pfam_family.isin(EPK | ELK) & m.cov_ok & m.len_ok]
with open(f"{K}/pkl_domains.fasta", "w") as fo:
    for r in cand.itertuples():
        fo.write(f">{r.target}\n{seqs[r.target][r.env_from - 1:r.env_to]}\n")
subprocess.run(f"diamond makedb --in {K}/refs/refs_kd.fasta -d {K}/refs/refs_kd --quiet", shell=True, check=True)
subprocess.run(f"diamond blastp -q {K}/pkl_domains.fasta -d {K}/refs/refs_kd -o {K}/pkl_vs_refs.tsv --ultra-sensitive "
               f"-k 25 -e 1e-3 --threads 8 --quiet --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore",
               shell=True, check=True)
bh = pd.read_csv(f"{K}/pkl_vs_refs.tsv", sep="\t", header=None,
                 names=["target", "ref", "pident", "alen", "qlen", "slen", "evalue", "bits"])
bh["ref_label"] = bh.ref.str.split("|").str[1]
bh["ref_gene"] = bh.ref.str.split("|").str[2]
bh = bh.sort_values(["target", "bits"], ascending=[True, False])
top = bh.groupby("target").head(1).set_index("target")
sec = bh[bh.ref_label != bh.target.map(top.ref_label)].groupby("target").head(1).set_index("target")
m["best_label"] = m.target.map(top.ref_label)
m["best_gene"] = m.target.map(top.ref_gene)
m["best_pident"] = m.target.map(top.pident)
m["best_bits"] = m.target.map(top.bits)
m["second_label"] = m.target.map(sec.ref_label)
m["second_bits"] = m.target.map(sec.bits)


def klass(r):
    f = r.pfam_family
    if f in EPK:
        if r.domain == "Eukaryota":
            if r.best_label == "TK" and r.group not in ("Metazoa", "Holozoa_other"):
                return "ePK_TKlike_nonholozoan"  # best hit human TK outside Holozoa (median 30 % identity): not a TK
            return f"ePK_{r.best_label}" if isinstance(r.best_label, str) and r.best_label not in ("RIO", "ABC1") \
                else "ePK_unassigned"
        return "prok_ePK" if r.domain in ("Bacteria", "Archaea") else "Asgard_ePK"
    if f in ELK:
        d = {"Eukaryota": "euk", "Bacteria": "bac", "Archaea": "arc", "Asgard": "asg"}[r.domain]
        return f"{f}_{d}"
    return f"otherPKL_{f}"


m["klass"] = m.apply(klass, axis=1)
cols = ["acc", "leaf", "taxid", "name", "domain", "group", "pfam_family", "pfam_acc", "klass", "best_label",
        "best_gene", "best_pident", "best_bits", "second_label", "second_bits", "dom_score", "n_pkl_domains",
        "arch_pkl", "tlen", "env_from", "env_to", "dom_len", "hmm_cov", "dom_plddt", "mean_plddt", "cov_ok",
        "plddt_ok", "len_ok", "pass", "target"]
m = m[cols].sort_values(["klass", "domain", "group", "name"])
m.to_csv(f"{OUT}/membership_all.tsv.gz", sep="\t", index=False)
s = (m.groupby(["pfam_family", "klass", "domain"]).agg(n_hits=("acc", "size"), n_cov_ok=("cov_ok", "sum"),
                                                       n_pass=("pass", "sum"), n_species=("leaf", "nunique"))
     .reset_index().sort_values(["pfam_family", "klass", "domain"]))
s.to_csv(f"{OUT}/membership_summary.tsv", sep="\t", index=False)
pd.set_option("display.width", 200); pd.set_option("display.max_rows", 400)
print(s[s.n_hits >= 5].to_string())
