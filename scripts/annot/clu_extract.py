#!/usr/bin/env python3
"""Step 1a: extract AFDB v3 cluster-file rows for panel proteins (scratch parquet files)."""
import duckdb, sys, time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; RAW=ROOT/"data/raw/afdb_clusters/v3"
OUT=Path(sys.argv[1]); OUT.mkdir(parents=True,exist_ok=True)
con=duckdb.connect(); con.execute("set threads=8; set memory_limit='60GB'")
con.execute(f"""create table panel as select acc from read_parquet('{ROOT}/data/processed/panel/protein_table.parquet')
 union select acc from read_parquet('{ROOT}/data/processed/panel/protein_table_prok.parquet')""")
print(con.execute("select count(*) from panel").fetchone(),flush=True)
t=time.time()
con.execute(f"""copy (select * from read_csv('{RAW}/5-allmembers-repId-entryId-cluFlag-taxId.tsv.gz',delim='\t',header=false,
  columns={{'repId':'VARCHAR','entryId':'VARCHAR','cluFlag':'INTEGER','taxId':'BIGINT'}}) where entryId in (select acc from panel))
  to '{OUT}/f5_panel.parquet' (format parquet)""")
print('f5',time.time()-t,flush=True)
