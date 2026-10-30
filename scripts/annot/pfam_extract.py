#!/usr/bin/env python3
"""Step 3a: Pfam 38.2 regions for panel proteins -> data/interim/annotations/pfam_regions_panel.parquet
Pfam-A.regions.tsv.gz columns: pfamseq_acc, seq_version, crc64, md5, pfamA_acc, seq_start, seq_end (envelope), ali_start, ali_end (HMM alignment)."""
import duckdb, time
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]; RAW = ROOT/"data/raw/pfam"
OUT = ROOT/"data/interim/annotations"; OUT.mkdir(parents=True, exist_ok=True)
con = duckdb.connect(); con.execute("set threads=4; set memory_limit='30GB'")
con.execute(f"""create table panel as select acc from read_parquet('{ROOT}/data/processed/panel/protein_table.parquet')
 union select acc from read_parquet('{ROOT}/data/processed/panel/protein_table_prok.parquet')""")
t = time.time()
con.execute(f"""copy (select pfamseq_acc acc, seq_version, md5 pfam_seq_md5, pfamA_acc pfam_acc, seq_start env_start, seq_end env_end, ali_start, ali_end
  from read_csv('{RAW}/Pfam-A.regions.tsv.gz',delim='\t',header=true,quote='',escape='',
   columns={{'pfamseq_acc':'VARCHAR','seq_version':'INTEGER','crc64':'VARCHAR','md5':'VARCHAR','pfamA_acc':'VARCHAR','seq_start':'INTEGER','seq_end':'INTEGER','ali_start':'INTEGER','ali_end':'INTEGER'}})
  where pfamseq_acc in (select acc from panel)) to '{OUT}/pfam_regions_panel.parquet' (format parquet)""")
print('pfam', time.time()-t)
