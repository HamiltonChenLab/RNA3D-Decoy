# RNA3D-Decoy analysis results

[Back to the repository README](../README.md)

This page reports the descriptive analyses of the RNA3D-Decoy prediction pool, including TM-score, RNA-only lDDT, RNA-only RMSD, and RNA–protein interface DockQ. It documents the analysis snapshot used for these results. The tables describe the observed decoys; they do not evaluate a model-quality assessment (QA) predictor.

## Contents

- [Analysis population and definitions](#analysis-population-and-definitions)
- [Dataset summary](#dataset-summary)
- [TM-score](#tm-score)
- [RNA-only lDDT](#rna-only-lddt)
- [RNA-only RMSD](#rna-only-rmsd)
- [DockQ for Level III](#dockq-for-level-iii)
- [Interpretation and limitations](#interpretation-and-limitations)
- [Provenance and numerical verification](#provenance-and-numerical-verification)

## Analysis population and definitions

The analysis includes **677 eligible targets**: 416 in Level I, 117 in Level II, and 144 in Level III. Level I contains single-chain RNA without metal ions; Level II contains single-chain RNA with metal ions; Level III contains RNA–protein complexes.

The target selection follows the archived *basic target information* and TM-score reference snapshot. Level I/II targets with RNA sequences shorter than 10 nucleotides were excluded: `1OSU`, `1RXA`, `2IXZ`, `2OJ7`, and `2OJ8`. `1QZA` is absent from that reference snapshot and is not in the analysis population. Level III has no sequence-length filtering. Thus, the analysis population differs from the 683 entries in the repository's split lists.

All available eligible decoys are retained, without downsampling, equal-count balancing, or missing-score imputation. AF3 contributes 31 structures per fully scored target, ESMFold2 contributes 32, and Protenix, Chai-1, and Boltz-2 contribute 30 each. The AF3 count includes its additional top-level model. ESMFold2 has no predictions for Level III target `8D97`.

| Metric / scoring profile | Level I | Level II | Level III | Total valid decoys |
| --- | --- | --- | --- | --- |
| TM-score | 63,648 | 17,901 | 22,000 | 103,549 |
| RNA-only lDDT, C4′ | 63,648 | 17,901 | 21,847 | 103,396 |
| RNA-only RMSD, C4′ | 63,648 | 17,901 | 20,623 | 102,172 |
| DockQ, RNA–protein interface | Not applicable | Not applicable | 22,000 | 22,000 |

**Scoring scope.** Level I/II TM-scores evaluate RNA using C3′ atoms. Level III TM-scores retain the existing RNA–protein complex scoring. The RNA-only lDDT and RMSD tables use RNA C4′ atoms at all three levels. The combined dataset summary below separately preserves the mixed-backbone Level III profile, which uses RNA C4′ and protein Cα atoms. Ion coordinates are excluded from these evaluations. Level II lDDT and RMSD use the rerun evaluations, superseding the older spreadsheet values.

**Statistical units.** “All decoys” pools every valid score within the stated method and level. “Target medians” first takes the median across a target's valid decoys for each method, then summarizes these medians across targets, giving each scored target equal weight. SD is the sample standard deviation (`ddof=1`), describing dispersion rather than standard error or a confidence interval. Higher TM-score, lDDT, and DockQ are better; lower RMSD is better.

## Dataset summary

**Table 1. Eligible targets, decoy counts, and all-decoy mean ± sample SD.**

| Level | No. eligible targets | Total decoys | Metric | AF3 | ESMFold2 | Protenix | Chai-1 | Boltz-2 |
|---|---|---|---|---|---|---|---|---|
| I | 416 | 63,648 | No. decoys | 12,896 | 13,312 | 12,480 | 12,480 | 12,480 |
|  |  |  | TM-score | 0.375 ± 0.152 | 0.335 ± 0.144 | 0.316 ± 0.152 | 0.289 ± 0.112 | 0.337 ± 0.147 |
|  |  |  | Global lDDT | 0.788 ± 0.115 | 0.761 ± 0.116 | 0.742 ± 0.110 | 0.719 ± 0.127 | 0.767 ± 0.112 |
|  |  |  | Global RMSD (Å) | 5.476 ± 7.458 | 6.225 ± 7.863 | 6.781 ± 7.793 | 8.157 ± 9.429 | 6.604 ± 8.319 |
| II | 117 | 17,901 | No. decoys | 3,627 | 3,744 | 3,510 | 3,510 | 3,510 |
|  |  |  | TM-score | 0.569 ± 0.207 | 0.431 ± 0.222 | 0.596 ± 0.202 | 0.283 ± 0.162 | 0.395 ± 0.238 |
|  |  |  | Global lDDT | 0.743 ± 0.120 | 0.677 ± 0.157 | 0.769 ± 0.102 | 0.572 ± 0.192 | 0.682 ± 0.173 |
|  |  |  | Global RMSD (Å) | 11.962 ± 6.805 | 23.253 ± 15.020 | 11.516 ± 6.422 | 30.186 ± 16.078 | 27.998 ± 17.898 |
| III | 144 | 22,000 | No. decoys | 4,464 | 4,576 | 4,320 | 4,320 | 4,320 |
|  |  |  | TM-score | 0.379 ± 0.386 | 0.294 ± 0.333 | 0.355 ± 0.380 | 0.281 ± 0.333 | 0.228 ± 0.272 |
|  |  |  | Global lDDT | 0.819 ± 0.118 | 0.752 ± 0.149 | 0.762 ± 0.146 | 0.700 ± 0.164 | 0.627 ± 0.198 |
|  |  |  | Global RMSD (Å) | 7.161 ± 8.991 | 11.594 ± 12.863 | 9.351 ± 8.029 | 14.031 ± 13.927 | 15.108 ± 14.743 |
| Total | 677 | 103,549 | No. decoys | 20,987 | 21,632 | 20,310 | 20,310 | 20,310 |

Levels I and II evaluate RNA. **Level III in this table evaluates the RNA–protein complex**: its lDDT and RMSD use C4′ + Cα. These Level III values must not be substituted for the RNA-only results in Tables 4 and 5. “No. eligible targets” counts targets in the selected population; ESMFold2 covers 143 of the 144 Level III targets.

## TM-score

**Table 2. Overall TM-score summaries using all eligible decoys and target medians.**

| Method | No. decoys | All decoys: mean ± SD | All decoys: median | Target medians: mean ± SD | Target medians: median |
|---|---|---|---|---|---|
| AF3 | 20,987 | 0.4091 ± 0.2422 | 0.3627 | 0.4058 ± 0.2376 | 0.3622 |
| ESMFold2 | 21,632 | 0.3429 ± 0.2157 | 0.3050 | 0.3406 ± 0.2123 | 0.3157 |
| Protenix | 20,310 | 0.3730 ± 0.2502 | 0.3049 | 0.3729 ± 0.2498 | 0.3021 |
| Chai-1 | 20,310 | 0.2862 ± 0.1893 | 0.2538 | 0.2859 ± 0.1882 | 0.2554 |
| Boltz-2 | 20,310 | 0.3238 ± 0.2045 | 0.2855 | 0.3239 ± 0.2038 | 0.2780 |

These overall summaries pool the RNA scores from Levels I/II with complex scores from Level III. AF3 has the highest pooled mean and median in this observed pool, while the level-specific comparisons show different method ordering, particularly in Level II, where Protenix has the highest mean. Overall summaries therefore do not establish a universal method ranking.

**Table 3. All-decoy TM-score threshold coverage (%).**

| Level | Statistic / threshold | AF3 | ESMFold2 | Protenix | Chai-1 | Boltz-2 |
|---|---|---|---|---|---|---|
| I | No. decoys | 12,896 | 13,312 | 12,480 | 12,480 | 12,480 |
|  | > 0.21 | 91.29% | 86.51% | 79.26% | 76.88% | 85.54% |
|  | > 0.45 | 22.92% | 15.32% | 11.17% | 6.98% | 15.82% |
|  | > 0.60 | 8.92% | 6.00% | 6.20% | 2.72% | 5.54% |
|  | > 0.80 | 3.03% | 2.22% | 3.00% | 0.42% | 2.90% |
| II | No. decoys | 3,627 | 3,744 | 3,510 | 3,510 | 3,510 |
|  | > 0.21 | 94.51% | 93.86% | 94.84% | 63.73% | 93.53% |
|  | > 0.45 | 78.55% | 29.46% | 82.82% | 12.51% | 23.50% |
|  | > 0.60 | 48.14% | 21.26% | 54.07% | 5.61% | 20.54% |
|  | > 0.80 | 10.28% | 10.95% | 13.05% | 1.77% | 13.28% |
| III | No. decoys | 4,464 | 4,576 | 4,320 | 4,320 | 4,320 |
|  | > 0.21 | 43.75% | 43.36% | 43.75% | 43.75% | 42.48% |
|  | > 0.45 | 42.67% | 23.27% | 40.12% | 23.59% | 13.61% |
|  | > 0.60 | 34.95% | 18.20% | 26.88% | 18.70% | 10.74% |
|  | > 0.80 | 23.86% | 15.76% | 22.89% | 15.32% | 9.44% |

Coverage is `100 × count(TM-score > threshold) / number of valid decoys` within each method and level. The comparison is strictly **greater than** the threshold, and the thresholds are cumulative, not mutually exclusive quality classes. The cutoffs are descriptive reference values, not independently validated success criteria for RNA structure prediction.

At a threshold of >0.60, AF3 has the largest fraction in Level I (8.92%) and Level III (34.95%), while Protenix leads Level II (54.07%). At >0.80, the Level II ordering changes: Boltz-2 reaches 13.28%, compared with 13.05% for Protenix. This endpoint dependence is visible only when the score distribution is considered beyond a single mean or median.

## RNA-only lDDT

**Table 4. RNA-only C4′ global lDDT by level and overall.**

| Level | Statistic | AF3 | ESMFold2 | Protenix | Chai-1 | Boltz-2 |
|---|---|---|---|---|---|---|
| I | No. valid decoys | 12,896 | 13,312 | 12,480 | 12,480 | 12,480 |
|  | No. scored targets | 416 | 416 | 416 | 416 | 416 |
|  | All decoys: mean ± SD | 0.7884 ± 0.1151 | 0.7615 ± 0.1161 | 0.7422 ± 0.1101 | 0.7190 ± 0.1266 | 0.7668 ± 0.1123 |
|  | All decoys: median | 0.8000 | 0.7715 | 0.7500 | 0.7392 | 0.7785 |
|  | Target medians: mean ± SD | 0.7886 ± 0.1142 | 0.7612 ± 0.1141 | 0.7423 ± 0.1099 | 0.7193 ± 0.1262 | 0.7672 ± 0.1112 |
|  | Target medians: median | 0.8002 | 0.7728 | 0.7500 | 0.7390 | 0.7785 |
| II | No. valid decoys | 3,627 | 3,744 | 3,510 | 3,510 | 3,510 |
|  | No. scored targets | 117 | 117 | 117 | 117 | 117 |
|  | All decoys: mean ± SD | 0.7430 ± 0.1197 | 0.6773 ± 0.1573 | 0.7691 ± 0.1021 | 0.5720 ± 0.1924 | 0.6824 ± 0.1727 |
|  | All decoys: median | 0.6839 | 0.6021 | 0.7296 | 0.4490 | 0.5745 |
|  | Target medians: mean ± SD | 0.7416 ± 0.1202 | 0.6788 ± 0.1539 | 0.7695 ± 0.1022 | 0.5724 ± 0.1924 | 0.6802 ± 0.1745 |
|  | Target medians: median | 0.6796 | 0.5927 | 0.7314 | 0.4476 | 0.5622 |
| III | No. valid decoys | 4,433 | 4,544 | 4,290 | 4,290 | 4,290 |
|  | No. scored targets | 143 | 142 | 143 | 143 | 143 |
|  | All decoys: mean ± SD | 0.8454 ± 0.1492 | 0.7536 ± 0.1603 | 0.7358 ± 0.1519 | 0.6568 ± 0.1840 | 0.7274 ± 0.1669 |
|  | All decoys: median | 0.8958 | 0.7679 | 0.7215 | 0.6500 | 0.7500 |
|  | Target medians: mean ± SD | 0.8466 ± 0.1498 | 0.7566 ± 0.1507 | 0.7378 ± 0.1480 | 0.6573 ± 0.1804 | 0.7310 ± 0.1630 |
|  | Target medians: median | 0.8994 | 0.7659 | 0.7227 | 0.6500 | 0.7500 |
| Overall | No. valid decoys | 20,956 | 21,600 | 20,280 | 20,280 | 20,280 |
|  | No. scored targets | 676 | 675 | 676 | 676 | 676 |
|  | All decoys: mean ± SD | 0.7926 ± 0.1279 | 0.7452 ± 0.1378 | 0.7455 ± 0.1194 | 0.6804 ± 0.1628 | 0.7438 ± 0.1409 |
|  | All decoys: median | 0.7993 | 0.7586 | 0.7387 | 0.7082 | 0.7640 |
|  | Target medians: mean ± SD | 0.7927 ± 0.1276 | 0.7459 ± 0.1334 | 0.7461 ± 0.1181 | 0.6808 ± 0.1614 | 0.7445 ± 0.1396 |
|  | Target medians: median | 0.8002 | 0.7580 | 0.7400 | 0.7067 | 0.7633 |

This profile contains **103,396 valid lDDT scores across 676 scored targets**. Target `1EMI` lacks C4′ scores for all methods, and ESMFold2 additionally lacks `8D97`. These are missing scores, not additional sequence-length exclusions. The score is `global_lDDT`, not `average_local_lDDT` or all-atom lDDT.

AF3 has the highest all-decoy mean lDDT overall and in Levels I and III. Protenix has the highest mean in Level II. Both all-decoy and target-median summaries are shown so that score pooling and equal-target summaries can be compared directly.

## RNA-only RMSD

**Table 5. RNA-only C4′ global RMSD (Å) by level and overall.**

| Level | Statistic | AF3 | ESMFold2 | Protenix | Chai-1 | Boltz-2 |
|---|---|---|---|---|---|---|
| I | No. valid decoys | 12,896 | 13,312 | 12,480 | 12,480 | 12,480 |
|  | No. scored targets | 416 | 416 | 416 | 416 | 416 |
|  | All decoys: mean ± SD | 5.4763 ± 7.4583 | 6.2250 ± 7.8632 | 6.7810 ± 7.7931 | 8.1573 ± 9.4286 | 6.6042 ± 8.3193 |
|  | All decoys: median | 2.5830 | 3.3325 | 3.9847 | 4.4363 | 3.4539 |
|  | Target medians: mean ± SD | 5.4470 ± 7.3942 | 6.1707 ± 7.7609 | 6.7623 ± 7.7881 | 8.1454 ± 9.3504 | 6.5956 ± 8.2940 |
|  | Target medians: median | 2.5414 | 3.2798 | 4.0167 | 4.5098 | 3.4419 |
| II | No. valid decoys | 3,627 | 3,744 | 3,510 | 3,510 | 3,510 |
|  | No. scored targets | 117 | 117 | 117 | 117 | 117 |
|  | All decoys: mean ± SD | 11.9624 ± 6.8046 | 23.2528 ± 15.0199 | 11.5163 ± 6.4223 | 30.1864 ± 16.0777 | 27.9981 ± 17.8977 |
|  | All decoys: median | 13.1394 | 25.8693 | 13.1694 | 40.5646 | 39.0179 |
|  | Target medians: mean ± SD | 11.8285 ± 6.0824 | 24.1732 ± 14.6021 | 11.4729 ± 6.3980 | 30.2852 ± 16.2034 | 28.0361 ± 17.6671 |
|  | Target medians: median | 13.2571 | 34.4461 | 13.0614 | 42.2992 | 40.4339 |
| III | No. valid decoys | 4,185 | 4,288 | 4,050 | 4,050 | 4,050 |
|  | No. scored targets | 135 | 134 | 135 | 135 | 135 |
|  | All decoys: mean ± SD | 3.1956 ± 4.2352 | 7.7411 ± 11.9494 | 5.1854 ± 4.9661 | 10.6963 ± 15.0058 | 9.7977 ± 15.5457 |
|  | All decoys: median | 1.3362 | 2.3058 | 3.4209 | 3.4913 | 2.5662 |
|  | Target medians: mean ± SD | 3.1578 ± 4.1286 | 7.8826 ± 12.0668 | 5.0524 ± 4.7568 | 10.6126 ± 14.8504 | 9.8603 ± 15.6327 |
|  | Target medians: median | 1.3873 | 2.3632 | 3.3506 | 3.5265 | 2.5721 |
| Overall | No. valid decoys | 20,708 | 21,344 | 20,040 | 20,040 | 20,040 |
|  | No. scored targets | 668 | 667 | 668 | 668 | 668 |
|  | All decoys: mean ± SD | 6.1514 ± 7.3713 | 9.5165 ± 12.1366 | 7.2879 ± 7.3628 | 12.5288 ± 14.6053 | 10.9967 ± 14.5244 |
|  | All decoys: median | 2.7925 | 3.5724 | 4.3408 | 4.8961 | 3.8873 |
|  | Target medians: mean ± SD | 6.1021 ± 7.1868 | 9.6725 ± 12.2062 | 7.2418 ± 7.3251 | 12.5218 ± 14.5738 | 11.0107 ± 14.4724 |
|  | Target medians: median | 2.7352 | 3.5056 | 4.3072 | 4.8555 | 3.8975 |

This profile contains **102,172 valid RMSD scores across 668 scored targets**. In Level III, 1,224 evaluations failed across eight targets because fewer than three matched atom pairs were available: `5MFU`, `5O3J`, `6N6A`, `6N6C`, `6N6D`, `6N6E`, `6N6G`, and `6N6H`. In addition, `1EMI` lacks scores for all methods, and ESMFold2 lacks `8D97`. Missing and failed evaluations are not assigned a score of zero.

AF3 has the lowest all-decoy mean and median RMSD overall. In Level II, Protenix has the lower mean (11.5163 Å versus 11.9624 Å for AF3), while AF3 has the slightly lower pooled median (13.1394 Å versus 13.1694 Å). These statistics summarize different features of the distribution.

**Table 6. Pairwise Spearman correlations (ρ) of target-level median RNA-only RMSD.**

| Predictor pair | Level I | Level II | Level III |
|---|---|---|---|
| AF3 – ESMFold2 | 0.931 | 0.641 | 0.795 |
| AF3 – Protenix | 0.867 | 0.857 | 0.731 |
| AF3 – Chai-1 | 0.818 | 0.588 | 0.729 |
| AF3 – Boltz-2 | 0.884 | 0.688 | 0.643 |
| ESMFold2 – Protenix | 0.889 | 0.710 | 0.785 |
| ESMFold2 – Chai-1 | 0.850 | 0.836 | 0.764 |
| ESMFold2 – Boltz-2 | 0.913 | 0.866 | 0.714 |
| Protenix – Chai-1 | 0.865 | 0.614 | 0.897 |
| Protenix – Boltz-2 | 0.913 | 0.709 | 0.873 |
| Chai-1 – Boltz-2 | 0.856 | 0.902 | 0.815 |

Each target's median uses all of its valid decoys. Each correlation uses targets shared by the predictor pair: 416 in Level I, 117 in Level II, and 135 in Level III, or 134 for Level III comparisons involving ESMFold2. Decoys from different predictors are not arbitrarily paired by file order.

Correlations are strongest and most consistent in Level I (ρ = 0.818–0.931), and more variable in Level II (0.588–0.902) and Level III (0.643–0.897). Higher correlation indicates more similar rankings of target error, not necessarily lower RMSD or greater structural similarity between the predictors' models.

## DockQ for Level III

**Table 7. All-decoy DockQ reference-quality distribution for RNA–protein complexes.**

| Method | No. decoys | Incorrect: D < 0.23 | Acceptable: 0.23 ≤ D < 0.49 | Medium: 0.49 ≤ D < 0.80 | High: D ≥ 0.80 |
|---|---|---|---|---|---|
| AF3 | 4,464 | 1,168 (26.16%) | 654 (14.65%) | 1,611 (36.09%) | 1,031 (23.10%) |
| ESMFold2 | 4,576 | 1,510 (33.00%) | 1,045 (22.84%) | 1,535 (33.54%) | 486 (10.62%) |
| Protenix | 4,320 | 2,193 (50.76%) | 883 (20.44%) | 642 (14.86%) | 602 (13.94%) |
| Chai-1 | 4,320 | 2,365 (54.75%) | 1,139 (26.37%) | 500 (11.57%) | 316 (7.31%) |
| Boltz-2 | 4,320 | 1,931 (44.70%) | 1,177 (27.25%) | 1,037 (24.00%) | 175 (4.05%) |
| Total | 22,000 | 9,167 (41.67%) | 4,898 (22.26%) | 5,325 (24.20%) | 2,610 (11.86%) |

Here, D denotes DockQ. These categories are mutually exclusive and follow the stated reference cutoffs. Values are counts and percentages of the available decoys for each method. Percentages may not sum to 100% because of rounding. The data cover Level III only; ESMFold2 has 143 scored targets and the other methods have 144. There is no Level III length filtering.

Across 22,000 Level III decoys, AF3 has the highest proportions of medium-or-higher reference quality (59.18%) and high reference quality (23.10%). ESMFold2 is second in medium-or-higher coverage (44.17%), while Protenix has a larger high-quality fraction than ESMFold2 (13.94% versus 10.62%). Overall, 41.67% fall below the 0.23 reference cutoff. These results describe the observed interface-score distribution under the stated cutoffs.

### DockQ score distributions

![DockQ distributions for all decoys, target medians, and target means](dockq_01_distributions.png)

**Figure D1. All-decoy and across-target DockQ distributions.** The left panel includes all 22,000 Level III decoys. The middle and right panels summarize the per-target median and mean DockQ for each method. Each target contributes one value per method to these two panels; ESMFold2 covers 143 targets and the other methods cover 144. The red horizontal reference is the pooled all-decoy median (0.3208). Higher DockQ indicates better agreement at the RNA–protein interface.

## Interpretation and limitations

- Decoys from the same target are related observations. Pooled distributions describe the generated pool; target-level comparisons use targets as the statistical unit.
- Method output counts differ. Best-observed scores, winner shares, and oracle comparisons reflect the actual observed generation budgets, not equal-budget comparisons.
- Level III complex TM-score and mixed-backbone metrics have a different evaluation scope from RNA-only metrics. Differences between levels also reflect differences in target composition and length.
- RMSD is calculated on matched atom pairs after a global superposition. It does not by itself penalize all unmatched residues; atom coverage should be considered when reusing these scores.
- Missingness is metric dependent. Comparisons requiring a shared population should explicitly align valid targets and decoys rather than assume that all four metrics have identical denominators.
- Score dispersion describes variation in reference similarity; it is not a direct measure of structural diversity. These analyses also do not establish QA-ranking performance or performance on unrelated future targets.

## Provenance and numerical verification

The page was prepared from the saved analysis CSV and the decoy-level inputs embedded in the existing lDDT, RMSD, and DockQ reports. Numerical verification recomputed the table summaries, TM-score threshold counts, DockQ category counts, and RMSD correlations from those saved inputs. It included **920 numerical checks** across the RNA-only and mixed-backbone summaries and the final threshold/correlation tables. The original structural scoring programs were not rerun for this page.

| Input snapshot | SHA-256 |
| --- | --- |
| Filtered TM-score decoy CSV | `6a090a99e920b006c9340d51d0fa1831f72fb515013581a7ee7d1f9cd42b60cc` |
| lDDT source CSV recorded by the report | `2931ed88e835cf623b893cc2b138a515a1420371f0c162bc0334c6688e1fb08d` |
| RMSD source CSV recorded by the report | `e30301710ad7cb81849b6c5bcaf41c564e17703f28e29945d5d06ad68c3f206e` |
| DockQ source CSV recorded by the report | `1c2e21ff8c72ce438e2255f85cddc1c15988cb30a7af65ed27022f058018fe4c` |

The source-CSV hashes for lDDT, RMSD, and DockQ identify the input snapshots recorded when the reports were produced. The present numerical verification used the preserved decoy records embedded in those reports. Missing software-version information is not inferred: the archived lDDT report does not record an OpenStructure version, and the DockQ producer configuration/version was not verified in the original analysis.

The full dataset is linked from the [repository README](../README.md#download-data).
