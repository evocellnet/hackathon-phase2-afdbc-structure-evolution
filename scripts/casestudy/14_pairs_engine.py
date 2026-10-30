"""X6 node 2, step 1b: the project's branch engine on pairwise distances (DESIGN 11.8: 'same branch engine').

Per superfamily (400 AFDB v4 models already on disk, data/interim/casestudy/pdb/{sf}/):
  1. store   strucevo.pairs.store_files.build_store_from_files (Foldseek createdb; AA/3Di/CA/pLDDT arrays)
  2. pairs   all 79,800 within-superfamily pairs, query = shorter protein; Foldseek structurealign
             --alignment-type 2 (3Di+AA = the project's 'tdiaa' alignment, D22) via strucevo.pairs.align.run_foldseek
  3. counts  strucevo.pairs.posthoc.pair_metrics (core = pLDDT >= 70 in both, AA_ORDER LUT for both alphabets)
  4. dists   aa_ml LG+G4(alpha 1.64), tdi_ml Q.3Di.AF+G4(alpha 1.60) (config/tuning.yaml), d_max 20; lddt_D with the
             tdiaa S0 table and SE model; SEs x kappa (tdiaa table, binned by aa_ml)
  5. engine  strucevo.stats.branch.fit_branches on the UFBoot>=95 reference topology (polytomies native),
             one 'HOG' = the superfamily, leaves = proteins (pseudo-taxids), se_method='model' (D2), nu = eta/2.
  Variant 'msa': the same engine on pairwise ML distances computed from the shared MAFFT L-INS-i columns
             (occupancy >= 50%, no pLDDT mask) -> isolates the effect of pairwise vs multiple alignment.
Outputs: data/interim/casestudy/pairs/{sf}_pairs.parquet; results/casestudy/engine_branches_{sf}_{variant}.tsv
Usage: python 14_pairs_engine.py <sf> [--threads 10]
"""
import argparse, json, os, sys, time
import multiprocessing as mp
import numpy as np, pandas as pd, polars as pl, yaml

R = "/scratch/spascare/agentic-ai/workspaces/phase2-afdbc-evolution"
sys.path.insert(0, f"{R}/src")
from strucevo.pairs.store_files import build_store_from_files
from strucevo.pairs.store import Store
from strucevo.pairs import posthoc
from strucevo.pairs.align import run_foldseek, cleanup_foldseek, FS_FIELDS
from strucevo.pairs.ml import ml_from_counts, lddt_D
from strucevo.distances import pair_counts, AA_ORDER
from strucevo.stats.trees import PanelTree
from strucevo.stats.branch import fit_branches

CS = f"{R}/data/interim/casestudy"
OUT = f"{R}/results/casestudy"
ap = argparse.ArgumentParser()
ap.add_argument("sf")
ap.add_argument("--threads", type=int, default=10)
a = ap.parse_args()
sf = a.sf
tun = yaml.safe_load(open(f"{R}/config/tuning.yaml"))
ALPHA = tun["alpha"]
sel = pd.read_csv(f"{OUT}/sample_members.tsv", sep="\t")
sel = sel[sel.superfamily == sf].reset_index(drop=True)
seqids = sel.seqid.tolist()
pid = {s: i + 1 for i, s in enumerate(seqids)}          # pseudo-taxid
t0 = time.time()

# ---------------------------------------------------------------- 1. store
sdir = f"{CS}/pairs/store_{sf}"
ents = [dict(path=f"{CS}/pdb/{sf}/{s}.pdb", name=s, acc=s, kind="af") for s in seqids]
info = build_store_from_files(ents, sdir, threads=4)
assert info["n_store"] == 400, info
st = Store(sdir)
idx = st.idx_of_name
key = st.key

# ---------------------------------------------------------------- 2-3. pairs + counts
pq = f"{CS}/pairs/{sf}_pairs.parquet"
if not os.path.exists(pq):
    L = {s: int(st.len[idx[s]]) for s in seqids}
    qa, ta = [], []
    for i in range(400):
        for j in range(i + 1, 400):
            a_, b_ = seqids[i], seqids[j]
            if (L[b_], b_) < (L[a_], a_):
                a_, b_ = b_, a_
            qa.append(a_); ta.append(b_)
    wd = f"{CS}/pairs/fs_{sf}"
    tsv, tm = run_foldseek(st.db, key[[idx[q] for q in qa]], key[[idx[t] for t in ta]], wd, "tdiaa", threads=a.threads)
    fs = pl.read_csv(tsv, separator="\t", has_header=False, new_columns=FS_FIELDS, quote_char=None,
                     schema_overrides={"lddtfull": pl.Utf8, "cigar": pl.Utf8, "query": pl.Utf8, "target": pl.Utf8})
    print(f"foldseek {tm} rows {fs.height} / {len(qa)}", flush=True)
    posthoc.STORE = st
    tasks = [(k, idx[r["query"]], idx[r["target"]], r["qstart"], r["tstart"], r["cigar"], r["lddtfull"],
              False, False) for k, r in enumerate(fs.iter_rows(named=True))]
    # check offsets convention against the pipeline (1-based qstart from convertalis)
    with mp.get_context("fork").Pool(a.threads) as pool:
        res = []
        for part in pool.imap_unordered(posthoc.worker_foldseek, [tasks[i:i + 200] for i in range(0, len(tasks), 200)]):
            res.extend(part)
    res.sort(key=lambda x: x[0])
    c_aa = np.stack([m["c_aa"] for _, m in res]); c_td = np.stack([m["c_tdi"] for _, m in res])
    rows = fs.select("query", "target", "alntmscore", "qtmscore", "ttmscore", "fident", "alnlen").to_pandas()
    rows["n_core"] = [m["n_core"] for _, m in res]
    rows["cov_short"] = [m["cov_short"] for _, m in res]
    rows["lddt_global"] = [m["lddt_global"] for _, m in res]
    rows["len_a"] = [int(st.len[idx[q]]) for q in rows["query"]]
    np.savez_compressed(f"{CS}/pairs/{sf}_counts.npz", c_aa=c_aa, c_tdi=c_td)
    cleanup_foldseek(wd)
    pl.from_pandas(rows).write_parquet(pq)
rows = pl.read_parquet(pq).to_pandas()
cz = np.load(f"{CS}/pairs/{sf}_counts.npz")
c_aa, c_td = cz["c_aa"], cz["c_tdi"]
print(f"pairs {len(rows)} ({time.time() - t0:.0f}s)", flush=True)


def kappa_scale(aa_ml, metric, aln="tdiaa"):
    kp = tun["kappa"][aln]
    e = np.asarray(kp["bin_edges"])
    b = np.clip(np.searchsorted(e, np.nan_to_num(aa_ml, nan=0.0), side="right") - 1, 0, len(kp[metric]) - 1)
    return np.asarray(kp[metric])[b]


def distances(c_aa, c_td, S=None, len_a=None, aln="tdiaa"):
    A = ml_from_counts(c_aa, "aa", alpha=ALPHA["aa"])
    T = ml_from_counts(c_td, "tdi", alpha=ALPHA["tdi"])
    d = pd.DataFrame({"aa_ml": A["ml"], "aa_ml_se": A["se"] * kappa_scale(A["ml"], "aa_ml", aln), "aa_sat": A["sat"],
                      "aa_p": A["p"], "n_sites": A["n_sites"],
                      "tdi_ml": T["ml"], "tdi_ml_se": T["se"] * kappa_scale(A["ml"], "tdi_ml", aln), "tdi_sat": T["sat"]})
    if S is not None:
        Ld = lddt_D(S, len_a, tun["s0"]["tdiaa"], tun["lddt_se_model"]["tdiaa"])
        d["lddt_D"] = Ld["lddt_D"]; d["lddt_D_se"] = Ld["lddt_D_se"]; d["lddt_sat"] = Ld["lddt_sat"]
    return d


D = distances(c_aa, c_td, rows.lddt_global.values, rows.len_a.values)
D["species_a"] = [pid[q] for q in rows["query"]]
D["species_b"] = [pid[t] for t in rows["target"]]
D["status"] = np.where((rows.n_core.values >= 30), "ok", "short")
D["hog"] = sf; D["split"] = "casestudy"
D = pd.concat([rows.reset_index(drop=True), D], axis=1)
D.to_parquet(f"{CS}/pairs/{sf}_dist_tdiaa.parquet")

# ---------------------------------------------------------------- MSA-column variant
def read_fa(p):
    d, k = {}, None
    for line in open(p):
        line = line.strip()
        if line.startswith(">"):
            k = line[1:].split()[0]; d[k] = []
        elif k:
            d[k].append(line)
    return {k: "".join(v) for k, v in d.items()}


aaM = read_fa(f"{CS}/aln/{sf}_linsi_aa.fasta"); tdM = read_fa(f"{CS}/aln/{sf}_linsi_3di.fasta")
qa = rows["query"].tolist(); ta = rows["target"].tolist()
cm_aa = pair_counts([aaM[q] for q in qa], [aaM[t] for t in ta], AA_ORDER)
cm_td = pair_counts([tdM[q] for q in qa], [tdM[t] for t in ta], AA_ORDER)
DM = distances(cm_aa, cm_td)
for c in ("species_a", "species_b", "hog", "split"):
    DM[c] = D[c].values
DM["status"] = np.where(DM.n_sites.values >= 30, "ok", "short")
DM.to_parquet(f"{CS}/pairs/{sf}_dist_msa.parquet")

# ---------------------------------------------------------------- 5. engine
import dendropy
tr = dendropy.Tree.get(path=f"{OUT}/{sf}_reference_topology_ufb95.nwk", schema="newick", preserve_underscores=True)
for lf in tr.leaf_node_iter():
    lf.taxon.label = str(pid[lf.taxon.label])
for n in tr.internal_nodes():
    n.label = None
panel = PanelTree.from_dendropy(tr)
inv = {v: k for k, v in pid.items()}
nu = {"tdi_ml": tun["eta"]["tdiaa"]["tdi_ml"] / 2, "lddt_D": tun["eta"]["tdiaa"]["lddt_D"] / 2}
summ = {}
for variant, dd, mets in (("tdiaa", D, ("aa_ml", "tdi_ml", "lddt_D")), ("msa", DM, ("aa_ml", "tdi_ml"))):
    bf = fit_branches(dd, panel, metrics=mets, nu=nu if variant == "tdiaa" else {}, aln=variant,
                      hog_species={sf: np.arange(1, 401)}, se_method="model", misfit=True)
    E = pd.DataFrame({k: v for k, v in bf.edges.items() if not k.startswith("_")})
    W = E.pivot_table(index="edge_id", columns="metric", values=["b", "b_se", "cov_aa"], aggfunc="first")
    W.columns = [f"{v}_{m}" for v, m in W.columns]
    W = W.reset_index()
    term = E.drop_duplicates("edge_id").set_index("edge_id")[["is_terminal", "share"]]
    W = W.join(term, on="edge_id")
    W["leaves"] = [";".join(sorted(inv[int(x)] for x in e.split(","))) for e in W.edge_id]
    W["n_leaves"] = [e.count(",") + 1 for e in W.edge_id]
    W = W.drop(columns="edge_id")
    W.to_csv(f"{OUT}/engine_branches_{sf}_{variant}.tsv", sep="\t", index=False)
    F = pd.DataFrame(bf.fits)
    summ[variant] = {"fits": F.drop(columns=[c for c in F.columns if c.startswith("_")]).to_dict("records"),
                     "n_pairs": int(len(dd)), "n_ok_aa": int(((dd.status == "ok") & ~dd.aa_sat & dd.aa_ml.notna()).sum()),
                     "n_ok_tdi": int(((dd.status == "ok") & ~dd.tdi_sat & dd.tdi_ml.notna()).sum()),
                     "n_hidden_para_flag": int(np.nansum(bf.pair_flags["hidden_para"]))}
summ["alpha"] = ALPHA; summ["nu"] = nu; summ["wall_s"] = round(time.time() - t0)
json.dump(summ, open(f"{OUT}/engine_summary_{sf}.json", "w"), indent=1, default=float)
print(json.dumps({k: (v if k != "tdiaa" and k != "msa" else {kk: vv for kk, vv in v.items() if kk != "fits"})
                  for k, v in summ.items()}, default=float, indent=1))
