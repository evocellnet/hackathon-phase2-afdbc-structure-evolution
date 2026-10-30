#!/usr/bin/env bash
# Pfam 38.2 per-UniProt regions + clans (EBI FTP, CC0) and TED v1 consensus domain summary (Zenodo 13908086, CC-BY-4.0)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
P=https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release
cd "$ROOT/data/raw/pfam"
for f in Pfam-A.regions.tsv.gz Pfam-A.clans.tsv.gz Pfam.version.gz md5_checksums pfamA.txt.gz; do curl -sS -L --retry 5 -C - -O "$P/$f"; done
md5sum Pfam-A.regions.tsv.gz Pfam-A.clans.tsv.gz pfamA.txt.gz > local_md5.txt
cd "$ROOT/data/raw/ted"
Z=https://zenodo.org/records/13908086/files
for f in ted_365m.domain_summary.cath.globularity.taxid.tsv.gz ted_214m_per_chain_segmentation.tsv.gz; do curl -sS -L --retry 5 -C - -o "$f" "$Z/$f?download=1"; done
md5sum *.tsv.gz > local_md5.txt
echo ALL_DONE
