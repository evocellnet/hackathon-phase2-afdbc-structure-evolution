# Deep-node alternative ages (E4 sensitivity): literature curator report (2026-10-30)

**Bottom line.** Strassert et al. 2021 gives point estimates and 95% HPDs for 10 of the 16 listed panel nodes (8 exact, 1 approximate, 1 inferred by taxonomy). Betts et al. 2018 gives only one panel-relevant interval, crown Eukaryota, because it reports CIs, not point estimates, and its eukaryote sampling is thin. Nothing is interpolated: nodes a source lacks are NA.

## Files

- `data/interim/timetree/deep_ages_alt.tsv`: all columns, including `timetree_median`. Stays under the TimeTree-ignored directory; never commit.
- `results/ages/deep_ages_literature.tsv`: literature columns only, plus the mapped Strassert node id, `mapping_type` and citations. No TimeTree values or text (checked by grep). Check that the "no age columns under `results/`" rule (ANALYSIS_PLAN V7) allows published-literature ages. `check_contract.py` was not found in the repo, so I could not test it. If V7 forbids any age column, move the file to `data/interim/` or rename the columns.

## Sources (both verified)

| | Citation | How I read it |
|---|---|---|
| S21 | Strassert JFH, Irisarri I, Williams TA, Burki F. *A molecular timescale for eukaryote evolution with implications for the origin of red algal-derived plastids.* Nat Commun 12:1879 (2021). doi:10.1038/s41467-021-22044-z, PMID 33767194, PMC7994803. Erratum: Nat Commun 12:3574, doi:10.1038/s41467-021-23847-w | Full text from PMC. Per-node numbers come from Suppl. Data 1 (`Molecular_clock_results.xlsx`, sheets "mean all nodes" and "HPD all nodes"), downloaded from the Springer ESM. |
| B18 | Betts HC, Puttick MN, Clark JW, Williams TA, Donoghue PCJ, Pisani D. *Integrated genomic and fossil evidence illuminates life's early evolution and eukaryote origin.* Nat Ecol Evol 2:1556-1562 (2018). doi:10.1038/s41559-018-0644-x, PMC6152910 | Full text from PMC and the SI PDF. Title wording: "early evolution and eukaryote origin" (not "history"). |

**S21 setting used.** This is the Fig. 3 analysis: MCMCTree, autocorrelated clock, 33 soft-bound **uniform** calibrations, Amorphea root, 320 genes, 136 OTUs.

- Posterior mean is in "mean all nodes" column I. The paper gives **means, not medians**, so `strassert2021_median` holds the posterior mean.
- Lo/hi are the 95% HPD from "HPD all nodes" block 1. The sheet labels the columns `sup95`/`inf95`, but the smaller value is the lower bound, so I took min and max.
- Each row's `source_detail` lists the range of means across the other analyses (clock model, calibration prior, root). That range is what the main text quotes, e.g. LECA 2386-1958 Ma, crown red algae 1984-1732 Ma.

**B18.** Only the main-text intervals exist ("joint 95% credibility intervals" across all clock and calibration analyses, Fig. 3 bars): crown-Eukaryota 1842-1210 Ma (Discussion; abstract "<1.84 Ga"). No point estimates are given in the text or SI, and the Bitbucket repo holds inputs only. `betts2018_median` is NA throughout.

## Mapping (panel node -> source node)

The mapping uses descendant panel species (names from `panel_tree_mapping.tsv`) against S21 OTU membership in the "mean all nodes" leaf lists.

| Panel node | Clade | TimeTree | S21 node: mean (95% HPD) | B18 | Mapping |
|---|---|---|---|---|---|
| N001 | Eukaryota crown | 1455 | n137: 2132 (1951-2335) | 1210-1842 | exact |
| N002 | Metamonada + Heterolobosea | 1220 | NA (no Metamonada) | NA | absent |
| N006 | Rhodophyta | 1156 | n213: 1821 (1725-1935) | NA | exact |
| N007 | Galdieria + Cyanidiococcus | 989 | NA (see below) | NA | absent |
| N008 | Viridiplantae | 974 | n193: 1112 (849-1397) | NA | exact |
| N009 | Embryophyta (ref. only, TT < 500) | 483 | n197: 534 (516-546) | NA | exact |
| N026 | Fungi | 875 | n259: 882 (663-1111) | NA | exact |
| N027 | Fungi minus Blastocladio/chytrids | 657 | NA | NA | absent |
| N033 | Ascomycota (Pezizo + Taphrino) | 528 | n262: 493 (366-676) | NA | exact |
| N044 | Holozoa | 948 | n243: 974 (805-1184) | NA | exact |
| N045 | Choanozoa | 892 | n244: 838 (735-969) | NA | exact |
| N046 | Metazoa (Trichoplax + Eumetazoa) | 637 | n246: 651 (634-672) | NA | approximate |
| N067 | Spiralia | 588 | NA | NA | absent |
| N068 | Lophotrochozoa | 550 | NA | NA | absent |
| N083 | Amoebozoa | 1206 | n265: 1026 (758-1652) | NA | taxonomy-inferred |
| N086 | SAR + Haptista + Cryptista | 942 | NA (non-monophyletic) | NA | absent |

**Notes on the judgement calls (all in `source_detail`).**

- **N007.** S21 has Cyanidioschyzon x Galdieria (n221: 1452, HPD 1292-1583) but not Cyanidiococcus. Cyanidioschyzonales branch earlier, so the panel MRCA is younger than n221. n221 is therefore only an upper bound.
- **N027, N067, N068.** S21 lacks Mucoromycota/Zoopagomycota, Brachiopoda, Platyhelminthes and Rotifera. Dikarya (n261: 646) and Capitella-Lottia-Octopus (n254: 564, HPD 555-574) are lower bounds only.
- **N046.** The panel node contains Trichoplax, which S21 lacks. Strictly it is Parahoxozoa, which lies between Metazoa (n245: 688, HPD 661-719) and Eumetazoa (n246). The row is flagged approximate.
- **N083.** Acanthamoeba is not an S21 OTU. I assume it falls in the Discosea-like branch opposite Dictyostelium, so the MRCA is the Amoebozoa crown. The HPD is very wide.
- **N086.** In S21 Cryptista is sister to Archaeplastida, and Haptista group with TSAR, so the panel clade is not monophyletic. Its MRCA in S21 is n139 (2028, HPD 1875-2195), which includes Archaeplastida. This is a topology conflict, not missing data. NA.
- **Fungi text number.** S21's text says fungi originate at "1759 to 1078 Ma". Suppl. Data 1 does not reproduce that for n259 (range of means 871-1229); the closest nodes are n258 Holomycota and n242 Opisthokonta. I used the spreadsheet values and did not use the text number.

## Important caveat for the E4 sensitivity (collapsed panel tree)

The panel tree collapsed Opisthokonta, Amorphea, Archaeplastida, SAR and others into the **N001 root polytomy**. Pairs whose MRCA is N001 in our tree (e.g. human vs yeast, Opisthokonta in S21) would receive the S21 root age (2132) if N001 is mapped as above, instead of their own node age (S21 Opisthokonta: 1359). To avoid this, assign S21 ages through the uncollapsed topology (`panel_tree_uncollapsed.nwk`), using the reference ages below. Otherwise N001 is only correct for pairs that truly span the eukaryote root.

Collapsed-out deep nodes, for reference (S21 mean, 95% HPD, node id; B18 interval where reported):

| Clade | S21 | B18 |
|---|---|---|
| Eukaryota root | n137: 2132 (1951-2335) | 1210-1842 |
| Amorphea | n238: 2087 (1911-2287) | |
| Obazoa (Opisthokonta + Apusozoa + Breviatea) | n240: 1790 (1527-2057) | |
| Opisthokonta | n242: 1359 (1051-1679) | 1125-1707 |
| Archaeplastida | n191: 1925 (1799-2063) | 1118-1667 |
| Cryptista + Archaeplastida | n190: 1963 (1829-2113) | |
| Diaphoretickes (CA + TSAR + Haptista) | n139: 2028 (1875-2195) | |
| TSAR / SAR (Alveolata + Stramenopiles + Rhizaria) | n143: 1741 (1593-1901) | 1115-1645 |
| Discoba (incl. Tsukubamonas) | n229: 1827 (1610-2054) | |

B18 plastid-origin interval (stem Archaeplastida) 1118-1774 Ma. S21 plastid origin (Archaeplastida) 2137-1807 mya is the range of means across analyses (n191). Holozoa, Metazoa and Bilateria are in the table above or in S21 n245 (688), n247 (616).

## Substantive observation

S21 is systematically **older** for the deepest nodes than TimeTree: root 2132 vs 1455, Rhodophyta 1821 vs 1156, and Archaeplastida/Viridiplantae. At the same time it is **close or younger** for several mid-depth nodes: Fungi 882 vs 875, Holozoa 974 vs 948, Choanozoa 838 vs 892, Amoebozoa 1026 vs 1206. B18 (1210-1842) overlaps TimeTree's root CI (1228-1683) and is the less extreme of the two. The E4 sensitivity will therefore mostly move the root and the plant/red-algal/discoban/SAR-side deep nodes, and will barely change the fungal and holozoan nodes.

## Limits

- No per-node B18 point estimates and no panel-level B18 nodes beyond the root. Figure 3 is bars only, and no data file was found. I did not digitise it.
- S21 reports posterior means (the column holds means), and its HPDs are partly prior-driven where calibrations are tight (Embryophyta, Metazoa, Eumetazoa).
- Taxonomic assignment of S21 OTU codes is from the abbreviated labels (e.g. `Ar_Colp`, `Am_Pessone`), not the supplementary OTU table, which I did not open (68 MB zip).
- No files committed. Temporary downloads are in the session scratchpad only.
