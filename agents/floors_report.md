# Noise floors (floor_random, floor_ident): report (analyst, 2026-10-30)

**Bottom line.**
- 100,000 `floor_random` and 17,768 `floor_ident` pairs ran through all three alignments, with rows out = rows in (117,768 per alignment, 0 `fail_align`; 1,787 s wall on 8 threads).
- Every random pair was re-verified after the build: 0 violations of species, cluster, Pfam, pLDDT or md5 (`data/interim/floors/random_pair_verification.json`).
- `results/floors/floors.json` and `figures/floors/floor_distributions.png` are written.
- Identical-sequence pairs are not identical structures. Under tdiaa the median TM is 0.950, and only 8.8% of pairs have TM ≥ 0.999. So η is clearly above 0.

## Pair construction (`scripts/floors/build_floor_pairs.py`; seed = SeedSequence(20261030).spawn(16)[3], i.e. stage `floor_random_pairs`)

**Candidate proteins**
- Eukaryote panel only.
- AFDB v4 model fetched, `seq_match` (md5 of the 2026 sequence = md5 of the v4 sequence), mean pLDDT ≥ 70, and 50 ≤ len ≤ 2700.
- 1,807,171 proteins.
- Length ≥ 50 is required because `n_core ≥ 50` is needed for status `ok`. Bin "<100" therefore means 50–99 aa.

**floor_random**
- Different species (`species_taxid`).
- Different AFDB v3 structural cluster. The cluster is `repId` in file 5 (`5-allmembers…`), which covers all 214M entries and so gives every protein an ID; singletons count as their own cluster.
- No shared Pfam unit, where unit = clan if the family has one, otherwise the family.
  - Pfam 38.2 `Pfam-A.regions` + `Pfam-A.clans`.
  - Proteins are mapped by accession and by sequence md5, since the Pfam sequence versions differ from AFDB v4.
  - Both proteins must carry ≥ 1 Pfam family. 78.9% of candidates do; by bin: 54% (<100), 74%, 79%, 82%, 82%, 83% (≥800).
- Partner length in [len_a, 2·len_a]: HOG-like length ratio (realised median ratio 1.35).
- Each protein is used ≤ 2 times, and there are no duplicate pairs.
- Design counts (<100 … ≥800): 10k, 20k, 20k, 20k, 18k, 12k = 100,000 pairs, every bin ≥ 2,000.
- 185,316 distinct proteins.
- Pair-level rejection counts during sampling were minimal: species 990, cluster 29, Pfam 2,882, usage cap 1,188.

**floor_ident**
- One random pair of distinct accessions per md5 group of identical v4 sequences, among candidates with pLDDT ≥ 70.
- **Available in the eukaryote panel: 64,038 md5 groups, 143,322 possible pairs.** Only 0.3% of the sampled pairs are within the same species, so these are almost all cross-species identical proteins.
- Sampled 3,000 groups per length bin (all 2,768 for ≥ 800) = 17,768 pairs.
- The ≥ 200 threshold is met by a factor of 700, so **no AFDB v4 strain extension was needed**.
- The prokaryote panel (pLDDT ≥ 70) has only 35 md5 groups / 139 pairs, so it was not used.

**Run:**
```
python scripts/floors/build_floor_pairs.py     # lists in data/interim/floors/ (pairs_floors.tsv, accessions.tsv, pair_build_stats.json)
python scripts/pairs/fetch_structures.py data/interim/floors/accessions.tsv data/interim/structures/floors --procs 4 --threads 64   # 216,733 entries
python scripts/pairs/run_pairs.py --pairs data/interim/floors/pairs_floors.tsv --foldcomp data/interim/structures/floors \
   --store data/interim/structures/floors_store --out data/interim/pairs_floors --n-shards 100 --threads 8
python scripts/floors/compute_floors.py
```
The output goes to `data/interim/pairs_floors/{seqonly,tmalign,tdiaa}/`, **not** `data/interim/pairs/`. The shards are named `shard_{k}` and would collide with the primary pair run; merge by `pair_set` when loading. ML distances use the placeholder α = 1 and should be recomputed with `recompute_ml.py` after tuning (`--s0 results/floors/floors.json` needs the table format below; see "Open items").

## Results

### S0: floor_random (bin = length of shorter protein; median / q95)
`S0` = all pairs with an alignment row. `S0_ok` = only status `ok`. This is the population on which `lddt_D` is evaluated, but it is conditioned on good coverage and core size, so the floor is higher. Both are in the JSON.

| bin | tdiaa LDDT | tdiaa TM | tmalign LDDT | tmalign TM | seqonly LDDT |
|---|---|---|---|---|---|
| <100 | 0.086 / 0.336 | 0.143 / 0.285 | 0.206 / 0.382 | 0.325 / 0.464 | 0.067 / 0.330 |
| 100–200 | 0.060 / 0.229 | 0.101 / 0.197 | 0.153 / 0.257 | 0.294 / 0.398 | 0.036 / 0.219 |
| 200–300 | 0.050 / 0.188 | 0.079 / 0.163 | 0.126 / 0.205 | 0.278 / 0.374 | 0.021 / 0.135 |
| 300–500 | 0.044 / 0.166 | 0.064 / 0.134 | 0.103 / 0.158 | 0.262 / 0.339 | 0.013 / 0.083 |
| 500–800 | 0.043 / 0.145 | 0.053 / 0.105 | 0.086 / 0.125 | 0.239 / 0.305 | 0.010 / 0.059 |
| ≥800 | 0.039 / 0.135 | 0.043 / 0.085 | 0.073 / 0.108 | 0.218 / 0.285 | 0.006 / 0.043 |

`S0_ok` (median / q95), LDDT:

| bin | tdiaa (n) | tmalign (n) | seqonly (n) |
|---|---|---|---|
| <100 | 0.276 / 0.533 (394) | 0.240 / 0.435 (820) | 0.264 / 0.523 (637) |
| 100–200 | 0.201 / 0.363 (1832) | 0.178 / 0.291 (6683) | 0.205 / 0.415 (1653) |
| 200–300 | 0.179 / 0.308 (1688) | 0.156 / 0.244 (5119) | 0.178 / 0.362 (768) |
| 300–500 | 0.168 / 0.266 (1561) | 0.140 / 0.205 (2267) | 0.162 / 0.364 (253) |
| 500–800 | 0.157 / 0.243 (1057) | 0.129 / 0.195 (409) | 0.162 / 0.424 (79) |
| ≥800 | 0.156 / 0.234 (647) | 0.144 / 0.220 (88) | 0.263 / 0.570 (38) |

**Sanity.**
- Random-pair TM of 0.2–0.3 holds for `tmalign` (medians 0.22–0.33).
- Under `tdiaa` TM is much lower (0.04–0.14). Foldseek 3Di+AA aligns only a small local part of random pairs, so TM normalised by len_a is small. The aligned fraction (`cov_short` median 0.17) is low too.
- LDDT is low everywhere (medians ≤ 0.2 under `tmalign`).
- Under random pairs `seqonly` has the lowest S0, because it is a global alignment with free end gaps.
- Statuses of random pairs:

| alignment | ok | fail_core | fail_cov |
|---|---|---|---|
| tdiaa | 7.2% | 62.9% | 29.9% |
| tmalign | 15.4% | 18.5% | 66.2% |
| seqonly | 3.4% | 91.4% | 5.1% |

  So `ok` random pairs are rare; this is the expected behaviour of the status rule.

### η, ν (floor_ident, status `ok`; tdiaa n = 17,004; tmalign n = 16,877; seqonly n = 17,004)

| alignment | η_tdi_ml (median) | ν_tdi | η for 1 − LDDT_global (raw) | ν | η_aa_ml | median TM | median LDDT_global |
|---|---|---|---|---|---|---|---|
| tdiaa | **0.0529** | 0.0264 | **0.0205** | 0.0102 | 0 | 0.950 | 0.980 |
| tmalign | 0.0547 | 0.0273 | 0.0421 | 0.0210 | 0 | 0.954 | 0.958 |
| seqonly | 0.0529 | 0.0264 | 0.0205 | 0.0102 | 0 | n/a | 0.980 |

- The tdiaa and seqonly η values are identical because both give the identity alignment for identical sequences (differences are only in the post-hoc statistics).
- The AFDB prediction noise is real: only 8.8% of identical pairs reach TM ≥ 0.999, 12.6% reach LDDT ≥ 0.999, and 3.1% have tdi_ml = 0.
- η_tdi_ml rises with protein length. By bin (<100 … ≥800), tdiaa: 0.037, 0.042, 0.046, 0.055, 0.063, 0.066 (tmalign: up to 0.075). η is length-dependent, which §1.4 treats as a single constant, so the PI may want a length-binned η. Values by bin are in `floors.json#<aln>.eta_by_bin`.
- The tmalign distribution has a heavy tail from misregistration (fail_cov 53, fail_core 838; the tail in tdi_ml and 1 − LDDT in the figure). This is consistent with D22. The medians are not affected.
- Identical-pair `fail_core` is 4.3% (low-pLDDT-core proteins).

## Caveats and open items for the PI
1. **S0 choice.** Which S0 should define `lddt_sat` and `lddt_D`: the unconditional table (`S0`) or the status-`ok` table (`S0_ok`)? For tdiaa the q95 differs by a factor of ~1.6–2 (e.g. 200–300: 0.188 vs 0.308). I put `S0` (unconditional) at the top of `floors.json`, with `S0_ok` as a second key. The `S0_ok` bins ≥ 500 for tmalign/seqonly have n < 100 and are not reliable.
2. **Not homology-proof.**
   - Pfam clan coverage: only 46% of Pfam-A families have a clan.
   - Unannotated proteins were excluded.
   - Remote homologues in unclustered, different families cannot be excluded. The floor therefore comes from "annotated" proteins, and a contamination by remote homologues would raise q95 slightly (conservative).
   - The shorter bin (<100) is the least represented (54% annotated).
3. The partner-length restriction (≤ 2 × len_a) and the usage cap (≤ 2) are my choices; the plan is silent on them.
4. `floors.json` contains η/ν and S0, but **not** ρ_m(bin), κ_m(bin) or the lddt_D SE model (§2.3, to be produced by the SE calibration), and `lddt_D`-based η needs S0 (here η is the raw 1 − LDDT_global, as instructed). The JSON layout is per alignment: `S0[bin].{LDDT_global,TM}.{n,median,q05,q95,q99}`, `eta_tdi_ml`, `nu_tdi_ml`, `eta_lddt_raw`, `nu_lddt_raw`, `eta_aa_ml`, plus status counts. `recompute_ml.py --s0` may need a small adapter to read this layout.
5. ML distances in these tables use α = 1 (placeholder). η_tdi_ml should be recomputed with `recompute_ml.py` once α_3Di is tuned.

## Files
- Code: `scripts/floors/build_floor_pairs.py`, `scripts/floors/compute_floors.py`
- Results: `results/floors/floors.json`, `figures/floors/floor_distributions.png`
- Data (not for git): `data/interim/floors/` (candidate/cluster/Pfam intermediates, pair lists, `pair_build_stats.json`, `random_pair_verification.json`), `data/interim/structures/floors{,_store}`, `data/interim/pairs_floors/`
- Log: `logs/floors_run_pairs.log`
