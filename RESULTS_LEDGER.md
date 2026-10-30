# RESULTS_LEDGER — one line per reported number
| date | quantity | value | script | commit | input md5 | reviewed by |
|---|---|---|---|---|---|---|
| 2026-10-30 | split.tsv md5 (20/80 HOG split; never regenerate) | 29da378aec763731b728b38fa262941f | scripts/hogs/build_hogs.py (seed 20261030) | 9a0ca11 | hogs.orthoxml 58b35ecafd20 | PI |
| 2026-10-30 | X4 (exploratory) n 1:1 ortholog pairs, residue-matched, both X-ray ≤3 Å cov ≥0.70, status ok (tdiaa): post-cutoff / pre-cutoff / any-date | 197 (160 HOGs) / 485 (358) / 795 (510) | scripts/h1e/analyse_h1e.py (seed 20261030/1010, B=2000) | 2a2e8ae | pair_meta.tsv 6186aa14b00d | pending |
| 2026-10-30 | X4 post tdiaa: mean d_3Di AF–AF / PDB–PDB (95% CI, HOG bootstrap) | 0.59 [0.52, 0.66] | scripts/h1e/analyse_h1e.py (seed 20261030/1010, B=2000) | 2a2e8ae | pair_meta.tsv 6186aa14b00d | pending |
| 2026-10-30 | X4 post tdiaa: mean Δ d_3Di = AF − PDB | −0.165 [−0.194, −0.139]; 93% of pairs Δ<0 | scripts/h1e/analyse_h1e.py (seed 20261030/1010, B=2000) | 2a2e8ae | pair_meta.tsv 6186aa14b00d | pending |
| 2026-10-30 | X4 post tdiaa: mean Δ (1 − global LDDT) | −0.043 [−0.050, −0.036]; ratio 0.70 [0.62, 0.76] | scripts/h1e/analyse_h1e.py (seed 20261030/1010, B=2000) | 2a2e8ae | pair_meta.tsv 6186aa14b00d | pending |
| 2026-10-30 | X4 post tdiaa: mean Δ (1 − TM short) | −0.026 [−0.035, −0.018]; ratio 0.69 [0.60, 0.78] | scripts/h1e/analyse_h1e.py (seed 20261030/1010, B=2000) | 2a2e8ae | pair_meta.tsv 6186aa14b00d | pending |
| 2026-10-30 | X4 post tdiaa: mean Δ core RMSD (Å) | −0.75 [−1.13, −0.44]; ratio 0.62 [0.51, 0.75] | scripts/h1e/analyse_h1e.py (seed 20261030/1010, B=2000) | 2a2e8ae | pair_meta.tsv 6186aa14b00d | pending |
| 2026-10-30 | X4 post tdiaa: slope of Δ d_3Di vs d_AA | +0.052 [+0.026, +0.073] (tune+none only: +0.014 [−0.034, +0.047]) | scripts/h1e/analyse_h1e.py (seed 20261030/1010, B=2000) | 2a2e8ae | pair_meta.tsv 6186aa14b00d | pending |
| 2026-10-30 | X4 post tdiaa: slope of Δ (1 − LDDT) vs d_AA | +0.012 [+0.005, +0.018] | scripts/h1e/analyse_h1e.py (seed 20261030/1010, B=2000) | 2a2e8ae | pair_meta.tsv 6186aa14b00d | pending |
| 2026-10-30 | X4 post tdiaa: slope of Δ (1 − TM) vs d_AA | +0.004 [−0.005, +0.013] | scripts/h1e/analyse_h1e.py (seed 20261030/1010, B=2000) | 2a2e8ae | pair_meta.tsv 6186aa14b00d | pending |
| 2026-10-30 | X4 post tdiaa: slope ratio d_3Di AF/PDB vs d_AA | 1.34 [1.15, 1.57] | scripts/h1e/analyse_h1e.py (seed 20261030/1010, B=2000) | 2a2e8ae | pair_meta.tsv 6186aa14b00d | pending |
| 2026-10-30 | X4 any-date tdiaa: mean d_3Di ratio AF/PDB; slope Δ d_3Di | 0.61 [0.58, 0.64]; +0.047 [+0.036, +0.057] | scripts/h1e/analyse_h1e.py (seed 20261030/1010, B=2000) | 2a2e8ae | pair_meta.tsv 6186aa14b00d | pending |
| 2026-10-30 | X4 floor (same protein AF vs own PDB, matched, tdiaa, post): median d_3Di / 1−LDDT / 1−TM / RMSD_core | 0.241 [0.232, 0.256] / 0.043 / 0.018 / 0.83 Å (n=357) | scripts/h1e/analyse_h1e.py (seed 20261030/1010, B=2000) | 2a2e8ae | pair_meta.tsv 6186aa14b00d | pending |
| 2026-10-30 | X4 floor pre-cutoff chains: median d_3Di / 1−TM | 0.235 [0.225, 0.243] / 0.015 (n=807) | scripts/h1e/analyse_h1e.py (seed 20261030/1010, B=2000) | 2a2e8ae | pair_meta.tsv 6186aa14b00d | pending |
| 2026-10-30 | X4 same 50 ortholog pairs, post − pre chains: mean ΔΔ d_3Di (tdiaa) | +0.050 [+0.015, +0.087] (AF bias smaller with post-cutoff structures) | scripts/h1e/analyse_h1e.py (seed 20261030/1010, B=2000) | 2a2e8ae | pair_meta.tsv 6186aa14b00d | pending |
| 2026-10-30 | X4 unmatched reference (node-1 style, tdiaa post): median Δ (1−TM) / Δ d_3Di | +0.009 / −0.139 (TM sign flips when residues are not matched) | scripts/h1e/analyse_h1e.py (seed 20261030/1010, B=2000) | 2a2e8ae | pair_meta.tsv 6186aa14b00d | pending |
