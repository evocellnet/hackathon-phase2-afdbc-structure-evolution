#!/usr/bin/env python3
"""Step 1: AFDB Clusters v3 annotation of panel proteins -> data/interim/annotations/afdb_cluster_annot.parquet

Semantics (verified empirically, see agents/annotations_report.md):
  file 5 (allmembers): repId = Foldseek-cluster representative of the entry's cluster (for AFDB50 members, cluFlag 1, the
     repId is the cluster rep of their AFDB50 representative); cluFlag 1 = AFDB50 member (non-representative), 2 = AFDB50
     representative that is a member of a non-singleton Foldseek cluster (present in file 1), 3 = fragment, 4 = singleton.
  A protein is in a NON-SINGLETON cluster iff its repId is a row of file 2 (2.30M clusters, each nMem >= 2).
  file 2: nMem = number of cluster members (AFDB50 reps, no fragments); LCAtaxId = cluster LCA.
Usage: afdb_cluster_annot.py <scratch_dir>
"""
import duckdb, sys, time
import polars as pl
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]; RAW = ROOT/"data/raw/afdb_clusters/v3"; TAX = ROOT/"data/raw/ncbi_taxonomy/taxdump"
OUT = ROOT/"data/interim/annotations"; OUT.mkdir(parents=True, exist_ok=True)
S = Path(sys.argv[1])
con = duckdb.connect(); con.execute("set threads=8; set memory_limit='60GB'")
con.execute(f"create table f5 as select * from read_parquet('{S}/f5_panel.parquet')")
con.execute(f"""create table f2 as select * from read_csv('{RAW}/2-repId_isDark_nMem_repLen_avgLen_repPlddt_avgPlddt_LCAtaxId.tsv.gz',delim='\t',header=false,
 columns={{'repId':'VARCHAR','isDark':'INT','nMem':'INT','repLen':'INT','avgLen':'DOUBLE','repPlddt':'DOUBLE','avgPlddt':'DOUBLE','LCA':'BIGINT'}})""")
# total number of AFDB entries (incl. AFDB50 non-rep members and fragments) assigned to each panel cluster rep, from full file 5
t = time.time()
con.execute(f"""create table tot as select repId, count(*) n_entries_total, sum((cluFlag=1)::int) n_afdb50_members, sum((cluFlag=3)::int) n_fragments
 from read_csv('{RAW}/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz',delim='\t',header=false,
 columns={{'repId':'VARCHAR','entryId':'VARCHAR','cluFlag':'INTEGER','taxId':'BIGINT'}})
 where repId in (select repId from f5 where repId in (select repId from f2)) group by 1""")
print('tot', time.time()-t, flush=True)
df = con.execute("""select f5.entryId acc, f5.repId afdb_rep_raw, f5.cluFlag afdb_cluflag,
  (f2.repId is not null) afdb_in_cluster, f2.repId afdb_cluster_rep, f2.nMem afdb_cluster_nmem, tot.n_entries_total afdb_cluster_nentries,
  f2.isDark::boolean afdb_dark, f2.repLen afdb_rep_len, f2.avgLen afdb_avg_len, f2.repPlddt afdb_rep_plddt, f2.avgPlddt afdb_avg_plddt, f2.LCA afdb_lca_taxid
  from f5 left join f2 on f2.repId=f5.repId left join tot on tot.repId=f5.repId""").to_arrow_table()
df = pl.from_arrow(df)
# taxonomy
def dmp(p): return [l.rstrip("\n").split("\t|\t") for l in open(p)]
nodes = {}
for r in dmp(TAX/"nodes.dmp"): nodes[int(r[0])] = (int(r[1]), r[2])
names = {}
for r in dmp(TAX/"names.dmp"):
    if r[3].startswith("scientific name"): names[int(r[0])] = r[1]
merged = {int(r[0]): int(r[1].split("\t")[0]) for r in dmp(TAX/"merged.dmp")}
def grp(t):
    t = merged.get(t, t); seen = 0
    while t in nodes and seen < 200:
        if t == 2759: return "Eukaryota"
        if t == 2: return "Bacteria"
        if t == 2157: return "Archaea"
        if t == 10239: return "Viruses"
        if t == 131567: return "cellular_organisms"
        p = nodes[t][0]
        if p == t: break
        t = p; seen += 1
    return "root_other"
lcas = [x for x in df["afdb_lca_taxid"].unique().to_list() if x is not None]
tab = pl.DataFrame({"afdb_lca_taxid": lcas,
  "afdb_lca_name": [names.get(merged.get(t, t)) for t in lcas],
  "afdb_lca_rank": [nodes.get(merged.get(t, t), (0, None))[1] for t in lcas],
  "afdb_lca_group": [grp(t) for t in lcas]}, schema_overrides={"afdb_lca_taxid": pl.Int64})
df = df.join(tab, on="afdb_lca_taxid", how="left")
df.write_parquet(OUT/"afdb_cluster_annot.parquet")
print(df.shape); print(df.group_by("afdb_cluflag","afdb_in_cluster").len().sort("afdb_cluflag","afdb_in_cluster"))
print(df.filter(pl.col("afdb_in_cluster")).group_by("afdb_lca_rank").len().sort("len",descending=True).head(8))
