"""X6 node 2, step 3: map selected branches to concrete structural differences.

For a branch (bipartition from results/casestudy/branch_ratios_{sf}_ml.tsv, oriented away from the outgroup leaf):
  clade  = the oriented side; sister = the other descendants of its parent in the outgroup-rooted ML tree
           (widened to the next ancestor until it has >= 5 sequences).
  reps   = highest-mean-pLDDT member of the modal subfamily of each side.
  1. US-align rep_clade vs rep_sister (TM-scores, RMSD, identity) and per-residue CA deviation after the US-align
     superposition (-m rotation matrix).
  2. Column view on the full MAFFT L-INS-i MSA: per column, total-variation distance between the clade and sister
     distributions of 3Di states (TV_3Di) and of AA states (TV_AA), occupancy >= 50 % on both sides.
     'Structure-shift columns' = TV_3Di >= 0.5. Consecutive shift columns (gap <= 2) form regions.
  3. Numbering: reference protein added with `mafft --add --keeplength` (actin: rabbit skeletal alpha-actin P68135,
     reported in mature numbering = UniProt - 2; tubulin: pig alpha-tubulin P02550, UniProt numbering).
     Known-region table (approximate literature boundaries; see REGIONS) + UniProt features of the reference
     (helix/strand/binding sites; UniProt REST, fetched once to data/interim/casestudy/refs/).
Outputs results/casestudy/structmap_{sf}_columns.tsv, structmap_{sf}_regions.tsv, structmap_summary.tsv (appended)
Usage: python 18_structural_map.py <sf> <branch_id> [<branch_id> ...]
"""
import sys, os, json, re, subprocess, urllib.request
import numpy as np, pandas as pd
from ete3 import Tree

R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
CS = f"{R}/data/interim/casestudy"
OUT = f"{R}/results/casestudy"
sf, bids = sys.argv[1], sys.argv[2:]
REF = {"actin": ("P68135", 2), "tubulin": ("P02550", 0)}[sf]
# approximate functional regions in the reference numbering (mature skeletal actin; pig alpha-tubulin)
REGIONS = {"actin": [("P-loop1 (S14-G15 phosphate loop)", 11, 16), ("D-loop (DNase I-binding loop)", 38, 52),
                     ("sensor loop (H73)", 70, 78), ("P-loop2 (G156-G158)", 154, 161), ("W-loop / target-binding cleft", 165, 172),
                     ("hydrophobic (H-)plug", 262, 274), ("C-terminus / subdomain 1 hydrophobic groove", 349, 375),
                     ("hinge / helix 137-144", 137, 144), ("SD2 base (60-69)", 60, 69)],
           "tubulin": [("T1 loop", 10, 15), ("H1-S2 loop", 33, 60), ("T2 loop", 69, 74), ("T3 loop", 98, 104),
                       ("T4 loop (GGGTGSG signature)", 140, 146), ("T5 loop", 170, 182), ("T6 loop", 204, 208),
                       ("T7 loop", 244, 251), ("H8 (catalytic E254 in alpha)", 252, 260), ("M-loop", 272, 288),
                       ("S9-S10 loop (alpha insertion; taxane site in beta)", 359, 372), ("H11-H12 / C-terminal helices", 385, 437),
                       ("C-terminal tail (E-hook)", 438, 451)]}[sf]


def rf(p):
    d, k = {}, None
    for line in open(p):
        line = line.strip()
        if line.startswith(">"):
            k = line[1:].split()[0]; d[k] = []
        elif k:
            d[k].append(line)
    return {k: "".join(v) for k, v in d.items()}


sel = pd.read_csv(f"{OUT}/sample_members.tsv", sep="\t")
sel = sel[sel.superfamily == sf].set_index("seqid")
leaves = sorted(sel.index)
bit = {l: 1 << i for i, l in enumerate(leaves)}
BR = pd.read_csv(f"{OUT}/branch_ratios_{sf}_ml.tsv", sep="\t").set_index("branch_id")
OG = json.load(open(f"{OUT}/branch_ratio_summary_{sf}_ml.json"))["outgroup_leaf_for_orientation"]
t = Tree(f"{CS}/iqtree/{sf}_linsi_aa_ml.treefile", format=1)
t.set_outgroup(t & OG)
for n in t.traverse("postorder"):
    n.m = bit[n.name] if n.is_leaf() else sum(c.m for c in n.children)
node_of = {n.m: n for n in t.traverse()}

# ---- reference numbering: add the reference to the full MSA
A = {k: v.upper() for k, v in rf(f"{CS}/aln/{sf}_linsi_aa.full.fasta").items()}
T3 = rf(f"{CS}/aln/{sf}_linsi_3di.full.fasta")
refs = rf(f"{CS}/refs/refs.fasta") if os.path.exists(f"{CS}/refs/refs.fasta") else {}
racc, off = REF
rkey = [k for k in refs if k.split("|")[0] == racc][0]
addp = f"{CS}/refs/{racc}_add.fasta"
if not os.path.exists(addp):
    open(f"{CS}/refs/{racc}.fa", "w").write(f">{racc}\n{refs[rkey]}\n")
    r = subprocess.run(["mafft", "--add", f"{CS}/refs/{racc}.fa", "--keeplength", "--thread", "4", "--quiet",
                        f"{CS}/aln/{sf}_linsi_aa.full.fasta"], capture_output=True, text=True, check=True)
    open(addp, "w").write(r.stdout)
rs = rf(addp)[racc].upper()
ncol = len(next(iter(A.values())))
assert len(rs) == ncol
refpos = np.full(ncol, -1)
p = 0
for j, c in enumerate(rs):
    if c != "-":
        p += 1
        refpos[j] = p - off   # mature numbering for skeletal actin
feat_p = f"{CS}/refs/{racc}_features.json"
if not os.path.exists(feat_p):
    try:
        js = urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{racc}.json", timeout=30).read()
        open(feat_p, "wb").write(js)
    except Exception as e:  # offline: no features
        json.dump({"features": []}, open(feat_p, "w"))
FE = []
for f in json.load(open(feat_p)).get("features", []):
    if f["type"] in ("Helix", "Beta strand", "Binding site", "Region", "Motif", "Site"):
        s, e = f["location"]["start"]["value"], f["location"]["end"]["value"]
        if s is None or e is None:
            continue
        lig = f.get("ligand", {}).get("name", "")
        FE.append((f["type"], s - off, e - off, (f.get("description") or lig)[:40]))
SS = {}
for ty, s, e, _ in FE:
    if ty in ("Helix", "Beta strand"):
        for q in range(s, e + 1):
            SS[q] = "H" if ty == "Helix" else "E"

keys = list(A)
MA = np.array([list(A[k]) for k in keys]); M3 = np.array([list(T3[k]) for k in keys])
kidx = {k: i for i, k in enumerate(keys)}


def tv(X):
    """per-column state frequency (non-gap) -> dict-free: return (freq matrix over states, occupancy)"""
    out = []
    return out


def coldist(rows1, rows2, M):
    d = np.full(ncol, np.nan); occ1 = (M[rows1] != "-").mean(0); occ2 = (M[rows2] != "-").mean(0)
    maj1 = np.empty(ncol, object); maj2 = np.empty(ncol, object)
    for j in range(ncol):
        a = M[rows1, j]; b = M[rows2, j]
        a = a[a != "-"]; b = b[b != "-"]
        if len(a) == 0 or len(b) == 0:
            continue
        sa, ca = np.unique(a, return_counts=True); sb, cb = np.unique(b, return_counts=True)
        fa = dict(zip(sa, ca / ca.sum())); fb = dict(zip(sb, cb / cb.sum()))
        d[j] = 0.5 * sum(abs(fa.get(s, 0) - fb.get(s, 0)) for s in set(fa) | set(fb))
        maj1[j] = sa[np.argmax(ca)]; maj2[j] = sb[np.argmax(cb)]
    return d, occ1, occ2, maj1, maj2


def usalign(p1, p2, mat):
    r = subprocess.run(["USalign", p1, p2, "-m", mat], capture_output=True, text=True, check=True).stdout
    tm = [float(x) for x in re.findall(r"TM-score= ([0-9.]+)", r)[:2]]
    rmsd = float(re.search(r"RMSD=\s*([0-9.]+)", r).group(1))
    ide = float(re.search(r"Seq_ID=n_identical/n_aligned=\s*([0-9.]+)", r).group(1))
    lines = [l for l in r.splitlines()]
    i = [k for k, l in enumerate(lines) if l.startswith('(":"')][0]
    s1, s2 = lines[i + 1], lines[i + 3]
    U = []; tvec = []
    for l in open(mat).read().splitlines()[2:5]:
        v = l.split()
        tvec.append(float(v[1])); U.append([float(x) for x in v[2:5]])
    return tm, rmsd, ide, s1, s2, np.array(U), np.array(tvec)


def ca(path):
    out = []
    for l in open(path):
        if l.startswith("ATOM") and l[12:16].strip() == "CA":
            out.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return np.array(out)


summ, regs, colrows = [], [], []
os.makedirs(f"{CS}/structmap", exist_ok=True)
for bid in bids:
    m = int(BR.at[bid, "clade_hex"], 16)
    node = node_of.get(m)
    if node is None:
        print("not a clade in the rooted tree", bid); continue
    anc = node.up
    sm = anc.m ^ node.m
    while bin(sm).count("1") < 5 and anc.up is not None:   # widen the sister side to >= 5 sequences
        anc = anc.up
        sm = anc.m ^ node.m
    C = [l for l in leaves if m & bit[l]]; S_ = [l for l in leaves if sm & bit[l]]
    def rep(ms):  # highest mean pLDDT among members of the side's modal subfamily
        s_ = sel.loc[ms]
        return s_[s_.subfamily == s_.subfamily.mode().iloc[0]].mean_plddt.idxmax()
    rc, rsis = rep(C), rep(S_)
    p1, p2 = f"{CS}/pdb/{sf}/{rc}.pdb", f"{CS}/pdb/{sf}/{rsis}.pdb"
    mat = f"{CS}/structmap/{bid}.mat"
    tm, rmsd, ide, s1, s2, U, tv_ = usalign(p1, p2, mat)
    X1 = ca(p1) @ U.T + tv_; X2 = ca(p2)
    # per-residue deviation of rep_clade residues (US-align alignment)
    dev1 = np.full(len(X1), np.nan)
    i1 = i2 = 0
    for a_, b_ in zip(s1, s2):
        if a_ != "-" and b_ != "-":
            dev1[i1] = np.linalg.norm(X1[i1] - X2[i2])
        i1 += a_ != "-"; i2 += b_ != "-"
    # MSA column of each rep_clade residue
    colof = np.flatnonzero(MA[kidx[rc]] != "-")
    devcol = np.full(ncol, np.nan); devcol[colof] = dev1
    rows1 = [kidx[k] for k in C]; rows2 = [kidx[k] for k in S_]
    d3, o1, o2, m31, m32 = coldist(rows1, rows2, M3)
    dA, _, _, mA1, mA2 = coldist(rows1, rows2, MA)
    ok = (o1 >= 0.5) & (o2 >= 0.5)
    shift = ok & (d3 >= 0.5)
    for j in np.flatnonzero(ok):
        colrows.append(dict(branch_id=bid, col=j, refpos=refpos[j], ss=SS.get(refpos[j], ""), occ_clade=o1[j], occ_sister=o2[j],
                            tv_3di=d3[j], tv_aa=dA[j], maj3di_clade=m31[j], maj3di_sister=m32[j], majaa_clade=mA1[j],
                            majaa_sister=mA2[j], ca_dev_reps=devcol[j], shift=bool(shift[j])))
    # regions
    js = np.flatnonzero(shift)
    groups = []
    for j in js:
        if groups and j - groups[-1][-1] <= 3:
            groups[-1].append(j)
        else:
            groups.append([j])
    for g in groups:
        rp = [refpos[j] for j in g if refpos[j] > 0]
        lo_, hi_ = (min(rp), max(rp)) if rp else (np.nan, np.nan)
        known = [n for n, a_, b_ in REGIONS if rp and not (hi_ < a_ or lo_ > b_)]
        fe = sorted({f"{ty}:{d}" if d else ty for ty, a_, b_, d in FE if rp and ty not in ("Helix", "Beta strand") and not (hi_ < a_ or lo_ > b_)})
        ss = "".join(SS.get(q, "-") for q in range(int(lo_), int(hi_) + 1)) if rp else ""
        regs.append(dict(branch_id=bid, clade_name=BR.at[bid, "clade_name"], ref=racc, ref_from=lo_, ref_to=hi_, n_cols=len(g),
                         mean_tv_3di=float(np.nanmean(d3[g])), mean_tv_aa=float(np.nanmean(dA[g])),
                         mean_ca_dev=float(np.nanmean(devcol[g])), ss_ref=ss, known_region=";".join(known),
                         uniprot_features=";".join(fe)[:200]))
    summ.append(dict(branch_id=bid, sf=sf, clade_name=BR.at[bid, "clade_name"], n_clade=len(C), n_sister=len(S_),
                     rep_clade=rc, rep_clade_taxon=sel.at[rc, "name"].split(" (")[0], rep_clade_sub=sel.at[rc, "subfamily"],
                     rep_sister=rsis, rep_sister_taxon=sel.at[rsis, "name"].split(" (")[0], rep_sister_sub=sel.at[rsis, "subfamily"],
                     tm_norm_clade=tm[0], tm_norm_sister=tm[1], rmsd=rmsd, seq_id=ide,
                     n_cols_ok=int(ok.sum()), n_shift_cols=int(shift.sum()),
                     frac_shift_cols=float(shift.sum() / max(ok.sum(), 1)), mean_tv_3di=float(np.nanmean(d3[ok])),
                     mean_tv_aa=float(np.nanmean(dA[ok])), n_cols_ca_dev_gt3=int(np.nansum(devcol[ok] > 3)),
                     rho_iq_linsi=BR.at[bid, "iq_linsi_rho"], b_aa_iq_linsi=BR.at[bid, "iq_linsi_b_aa"]))
    print(summ[-1])
pd.DataFrame(colrows).to_csv(f"{OUT}/structmap_{sf}_columns.tsv", sep="\t", index=False)
pd.DataFrame(regs).to_csv(f"{OUT}/structmap_{sf}_regions.tsv", sep="\t", index=False)
S = pd.DataFrame(summ)
sp = f"{OUT}/structmap_summary.tsv"
if os.path.exists(sp):
    old = pd.read_csv(sp, sep="\t"); old = old[old.sf != sf]; S = pd.concat([old, S])
S.to_csv(sp, sep="\t", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 60)
print(pd.DataFrame(regs).sort_values(["branch_id", "mean_tv_3di"], ascending=[True, False]).groupby("branch_id").head(8).to_string())

# per known functional region: column statistics vs the whole-core average (descriptive)
Cdf = pd.DataFrame(colrows)
kr = []
for bid, d in Cdf.groupby("branch_id"):
    base = dict(tv3=d.tv_3di.mean(), tva=d.tv_aa.mean(), dev=d.ca_dev_reps.mean(), sh=d["shift"].mean())
    for nm, a_, b_ in REGIONS:
        r = d[(d.refpos >= a_) & (d.refpos <= b_)]
        if len(r) == 0:
            continue
        kr.append(dict(branch_id=bid, region=nm, ref_from=a_, ref_to=b_, n_cols=len(r), frac_shift=r["shift"].mean(),
                       frac_shift_core=base["sh"], mean_tv_3di=r.tv_3di.mean(), mean_tv_3di_core=base["tv3"],
                       mean_tv_aa=r.tv_aa.mean(), mean_tv_aa_core=base["tva"], mean_ca_dev=r.ca_dev_reps.mean(),
                       mean_ca_dev_core=base["dev"]))
KR = pd.DataFrame(kr)
KR.to_csv(f"{OUT}/structmap_{sf}_known_regions.tsv", sep="\t", index=False)
print(KR.round(2).to_string())
