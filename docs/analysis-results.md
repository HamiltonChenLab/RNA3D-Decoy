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

![T1: All-decoy TM-score threshold exceedance by method and level](TMScore-classification.png)

**Figure T1.** All-decoy TM-score threshold coverage corresponding to Table 3. Each cell is the percentage of valid decoys strictly above 0.21, 0.45, 0.60, or 0.80 within one method and level, rounded to one decimal place. These are cumulative threshold exceedances, not mutually exclusive quality categories.

### All-decoy TM-score figures

All figures use the analysis population above. Level I/II TM-scores evaluate RNA; Level III TM-scores evaluate RNA–protein complexes. Click an image to inspect its full resolution.

![T2: fig 02d method all decoy violin ](fig_02d_method_all_decoy_violin.png)

**Figure T2.** All-decoy TM-score distributions by method. Every eligible observed decoy is included; the overlaid boxes summarize the pooled distributions.

![T3: fig raw 01 method level violin ](fig_raw_01_method_level_violin.png)

**Figure T3.** All-decoy TM-score distributions by method and level. Every eligible observed decoy is retained; the inner boxes show the median and interquartile range within each method and level. Level I/II scores evaluate RNA, whereas Level III scores evaluate RNA–protein complexes.

![T4: fig raw 02 method level ecdf ](fig_raw_02_method_level_ecdf.png)

**Figure T4.** Exact all-decoy empirical cumulative distributions by method within each level. Curves show the fraction of observed decoys with TM-score at or below each value.

![T5: fig all 06 sequence length vs tm score ](fig_all_06_sequence_length_vs_tm_score.png)

**Figure T5.** All-decoy TM-score versus sequence length, separated by method and level. Level III length includes RNA and protein in the evaluated complex; its scale is not RNA-only length.

![T6: fig raw 04 target score landscape by method ](fig_raw_04_target_score_landscape_by_method.png)

**Figure T6.** TM-score distributions across all targets for each method. All observed decoys are retained, with target-level summary intervals as shown in the figure.

![T7: fig raw 05 target method threshold heatmaps ](fig_raw_05_target_method_threshold_heatmaps.png)

**Figure T7.** Within-target, within-method fractions of decoys with TM-score strictly above 0.20, 0.40, 0.60, and 0.80. These exploratory cutoffs differ from Table 3, which uses 0.21, 0.45, 0.60, and 0.80; missing cells are not zero.

<details>
<summary>Additional TM-score distributions, target summaries, and method complementarity (18 figures)</summary>

These analyses complement the all-decoy views. Bootstrap error bars are 95% pointwise percentile intervals from 10,000 target-level resamples (seed 20260901). They do not represent simultaneous confidence intervals or equal generation budgets.

![T8: fig 01a overall tm histogram ](fig_01a_overall_tm_histogram.png)

**Figure T8.** Pooled TM-score histogram for all 103,549 decoys, using bins of width 0.02. Each observed decoy contributes equally.

![T9: fig 01b method ecdf ](fig_01b_method_ecdf.png)

**Figure T9.** Exact pooled empirical cumulative distributions for each method's observed decoys. Unequal generation counts are retained.

![T10: fig 01c level ecdf ](fig_01c_level_ecdf.png)

**Figure T10.** Exact pooled empirical cumulative distributions by level. Levels contain different targets, and Level III evaluates RNA–protein complexes.

![T11: fig 02a method target median ](fig_02a_method_target_median.png)

**Figure T11.** Across-target distributions of target–method median TM-scores. Each observation represents one target–method median.

![T12: fig 02b method target best ](fig_02b_method_target_best.png)

**Figure T12.** Across-target distributions of the best observed TM-score in each target–method cell. These maxima depend on the actual numbers of generated decoys.

![T13: fig 02c method quality breadth ](fig_02c_method_quality_breadth.png)

**Figure T13.** Within-target TM-score interquartile ranges by method. Score dispersion describes quality breadth, not three-dimensional structural diversity.

![T14: fig 03a target difficulty ](fig_03a_target_difficulty.png)

**Figure T14.** Target-level pooled median TM-scores, combining the available methods for each target. Difficulty is relative to this observed method pool.

![T15: fig 03b target median vs iqr ](fig_03b_target_median_vs_iqr.png)

**Figure T15.** Target-level median TM-score versus within-target interquartile range, grouped by level. Each point is one target.

![T16: fig 03c target quality entropy ](fig_03c_target_quality_entropy.png)

**Figure T16.** Normalized entropy across the descriptive TM-score bins [0, 0.2), [0.2, 0.4), [0.4, 0.6), [0.6, 0.8), and [0.8, 1]. These bins are not validated folding-success categories.

![T17: fig 03d sorted target intervals ](fig_03d_sorted_target_intervals.png)

**Figure T17.** All targets sorted by their pooled median TM-score, with Q10–Q90 and Q25–Q75 intervals. No targets are subsampled.

![T18: fig 04a target method heatmap ](fig_04a_target_method_heatmap.png)

**Figure T18.** Target–method median TM-score heatmap. Rows represent targets and columns represent methods; missing observations remain marked as missing.

![T19: fig 04b method correlation ](fig_04b_method_correlation.png)

**Figure T19.** Pairwise Spearman correlations of target–method median TM-scores, using jointly observed targets for each pair. Score correlation does not establish structural redundancy.

![T20: fig 05a exact winner share ](fig_05a_exact_winner_share.png)

**Figure T20.** Shares of targets for which each method supplies the highest observed TM-score. Ties split credit equally, the denominator is all targets, and error bars are target-bootstrap intervals.

![T21: fig 05b near best coverage ](fig_05b_near_best_coverage.png)

**Figure T21.** Fractions of available targets for which a method's best observed score is within 0.01, 0.02, or 0.05 of the pooled best observed score.

![T22: fig 05c leave one out oracle drop ](fig_05c_leave_one_out_oracle_drop.png)

**Figure T22.** Mean loss in the pooled best observed TM-score after removing each method, with target-bootstrap intervals. Comparisons require the method and at least one alternative to be available.

![T23: fig 06a method level median ](fig_06a_method_level_median.png)

**Figure T23.** Medians across target–method median TM-scores within each method and level, with bootstrap intervals obtained by resampling targets within levels.

![T24: fig 06b method level best ](fig_06b_method_level_best.png)

**Figure T24.** Medians across target–method best observed TM-scores within each method and level, with within-level target-bootstrap intervals.

![T25: fig 06c method level quality breadth ](fig_06c_method_level_quality_breadth.png)

**Figure T25.** Medians across within-target TM-score interquartile ranges for each method and level, with within-level target-bootstrap intervals.

</details>

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

### RNA-only lDDT figures

These figures use C4′ atoms in all three levels. Higher lDDT indicates better local distance agreement. Counts and missing scores follow Table 4.

![L1: All decoys and target medians: violin distributions ](lddt_rna_only_01_distributions.png)

**Figure L1.** All-decoy and across-target-median lDDT distributions, shown overall and by level. Target medians are calculated within each target–method cell before being compared across targets.

![L2: All decoys and target medians: exact ECDFs ](lddt_rna_only_02_ecdf.png)

**Figure L2.** Exact empirical cumulative distributions of all-decoy lDDT and target-median lDDT. All available observations are retained.

![L3: All decoys and target medians: strict threshold coverage ](lddt_rna_only_03_threshold_coverage.png)

**Figure L3.** Fractions above the descriptive lDDT cutoffs shown in the figure. All-decoy and target-median coverage use different observation units; the comparisons are strict (>).

![L4: Target×method mean, SD, median and best ](lddt_rna_only_04_cell_summaries.png)

**Figure L4.** Target–method mean, sample SD, median, and best observed lDDT. Best observed means the maximum available score.

![L5: Target×method raw threshold fractions ](lddt_rna_only_05_cell_thresholds.png)

**Figure L5.** Within-target, within-method fractions of decoys strictly above the displayed lDDT cutoffs. Missing cells remain missing.

![L6: Target landscape and within-target quality breadth ](lddt_rna_only_06_target_landscape.png)

**Figure L6.** Target-level lDDT distributions and within-target quality breadth across the observed prediction pool.

![L7: Observed-pool method complementarity ](lddt_rna_only_07_complementarity.png)

**Figure L7.** Complementarity of methods in the observed prediction pool. Best-score and oracle comparisons reflect the actual generation counts for each method.

![L8: Across-target method correlations ](lddt_rna_only_08_correlations.png)

**Figure L8.** Across-target correlations between methods' score summaries. These describe shared score patterns and do not establish structural equivalence.

<details>
<summary>Mixed-backbone lDDT companion figures (8 figures)</summary>

This companion profile uses C4′ in Levels I/II and C4′ + Cα in Level III. Its Level III scores evaluate complexes and must be kept separate from the RNA-only results in Table 4. The profile corresponds to the scoring scope of the dataset summary in Table 1.

![L9: All decoys and target medians: violin distributions ](lddt_sheet_reference_01_distributions.png)

**Figure L9.** All-decoy and across-target-median lDDT distributions, shown overall and by level. Target medians are calculated within each target–method cell before being compared across targets. This figure uses the mixed-backbone companion profile.

![L10: All decoys and target medians: exact ECDFs ](lddt_sheet_reference_02_ecdf.png)

**Figure L10.** Exact empirical cumulative distributions of all-decoy lDDT and target-median lDDT. All available observations are retained. This figure uses the mixed-backbone companion profile.

![L11: All decoys and target medians: strict threshold coverage ](lddt_sheet_reference_03_threshold_coverage.png)

**Figure L11.** Fractions above the descriptive lDDT cutoffs shown in the figure. All-decoy and target-median coverage use different observation units; the comparisons are strict (>). This figure uses the mixed-backbone companion profile.

![L12: Target×method mean, SD, median and best ](lddt_sheet_reference_04_cell_summaries.png)

**Figure L12.** Target–method mean, sample SD, median, and best observed lDDT. Best observed means the maximum available score. This figure uses the mixed-backbone companion profile.

![L13: Target×method raw threshold fractions ](lddt_sheet_reference_05_cell_thresholds.png)

**Figure L13.** Within-target, within-method fractions of decoys strictly above the displayed lDDT cutoffs. Missing cells remain missing. This figure uses the mixed-backbone companion profile.

![L14: Target landscape and within-target quality breadth ](lddt_sheet_reference_06_target_landscape.png)

**Figure L14.** Target-level lDDT distributions and within-target quality breadth across the observed prediction pool. This figure uses the mixed-backbone companion profile.

![L15: Observed-pool method complementarity ](lddt_sheet_reference_07_complementarity.png)

**Figure L15.** Complementarity of methods in the observed prediction pool. Best-score and oracle comparisons reflect the actual generation counts for each method. This figure uses the mixed-backbone companion profile.

![L16: Across-target method correlations ](lddt_sheet_reference_08_correlations.png)

**Figure L16.** Across-target correlations between methods' score summaries. These describe shared score patterns and do not establish structural equivalence. This figure uses the mixed-backbone companion profile.

</details>

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

![R1: Target-median RMSD correlations and head-to-head comparisons across levels](fig_rmsd_11_median_head_to_head_vertical.png)

**Figure R1.** Pairwise comparisons of target-level median RNA-only RMSD, complementing Table 6. Upper-triangle cells show Spearman correlations and pairwise target counts; lower-triangle panels compare the same targets between two methods, with dashed identity lines. Marginal histograms show the target-median distributions. The displayed RMSD axes span 0–50 Å. Level I and II pairs include 416 and 117 targets; Level III pairs include 135 targets, or 134 when ESMFold2 is involved.

### RNA-only RMSD figures

These figures use C4′ atoms in all three levels and include every successfully evaluated eligible decoy. RMSD is reported in Å; lower is better. Long tails are retained. The target-median correlations are presented in Table 6 and Figure R1.

![R2: All-decoy and target-median violin distributions ](rmsd_rna_only_01_distributions.png)

**Figure R2.** All-decoy and across-target-median RMSD distributions, overall and by level. Lower RMSD indicates closer agreement after the stated structural superposition.

![R3: Mean ± sample SD and median ](rmsd_rna_only_02_mean_sd.png)

**Figure R3.** Mean, sample SD, and median RMSD for the observation units shown. SD describes dispersion; it is not a standard error or confidence interval.

![R4: Exact cumulative distributions ](rmsd_rna_only_03_ecdf.png)

**Figure R4.** Exact empirical cumulative distributions of RMSD. Every valid observation contributes, including the long upper tails.

![R5: Strict low-RMSD threshold coverage ](rmsd_rna_only_04_thresholds.png)

**Figure R5.** Fractions with RMSD strictly below 2, 5, 10, or 20 Å. These are descriptive cutoffs rather than validated success criteria.

![R6: Target×method mean, SD, median and minimum ](rmsd_rna_only_05_cell_summaries.png)

**Figure R6.** Target–method mean, sample SD, median, and minimum observed RMSD. The minimum is the best observed score for this metric.

![R7: Target×method raw threshold fractions ](rmsd_rna_only_06_cell_thresholds.png)

**Figure R7.** Within-target, within-method fractions of valid decoys below the displayed RMSD cutoffs. Failed or absent scores are not replaced by zero.

![R8: Target landscape and error dispersion ](rmsd_rna_only_07_landscape.png)

**Figure R8.** Target-level RMSD distributions and within-target error dispersion. Any symlog axis is a display transform; statistics retain the original Å values.

![R9: Observed-pool method complementarity ](rmsd_rna_only_08_complementarity.png)

**Figure R9.** Observed-pool method complementarity using minimum RMSD as the best score. Generation counts and the available target–method cells determine the observed oracle.

<details>
<summary>Mixed-backbone RMSD companion figures (8 figures)</summary>

This companion profile uses C4′ in Levels I/II and a shared C4′ + Cα superposition for Level III RNA–protein complexes. Its Level III values differ in scope from RNA-only RMSD and correspond to the dataset-summary convention in Table 1.

![R10: All-decoy and target-median violin distributions ](rmsd_mixed_backbone_01_distributions.png)

**Figure R10.** All-decoy and across-target-median RMSD distributions, overall and by level. Lower RMSD indicates closer agreement after the stated structural superposition. This figure uses the mixed-backbone companion profile.

![R11: Mean ± sample SD and median ](rmsd_mixed_backbone_02_mean_sd.png)

**Figure R11.** Mean, sample SD, and median RMSD for the observation units shown. SD describes dispersion; it is not a standard error or confidence interval. This figure uses the mixed-backbone companion profile.

![R12: Exact cumulative distributions ](rmsd_mixed_backbone_03_ecdf.png)

**Figure R12.** Exact empirical cumulative distributions of RMSD. Every valid observation contributes, including the long upper tails. This figure uses the mixed-backbone companion profile.

![R13: Strict low-RMSD threshold coverage ](rmsd_mixed_backbone_04_thresholds.png)

**Figure R13.** Fractions with RMSD strictly below 2, 5, 10, or 20 Å. These are descriptive cutoffs rather than validated success criteria. This figure uses the mixed-backbone companion profile.

![R14: Target×method mean, SD, median and minimum ](rmsd_mixed_backbone_05_cell_summaries.png)

**Figure R14.** Target–method mean, sample SD, median, and minimum observed RMSD. The minimum is the best observed score for this metric. This figure uses the mixed-backbone companion profile.

![R15: Target×method raw threshold fractions ](rmsd_mixed_backbone_06_cell_thresholds.png)

**Figure R15.** Within-target, within-method fractions of valid decoys below the displayed RMSD cutoffs. Failed or absent scores are not replaced by zero. This figure uses the mixed-backbone companion profile.

![R16: Target landscape and error dispersion ](rmsd_mixed_backbone_07_landscape.png)

**Figure R16.** Target-level RMSD distributions and within-target error dispersion. Any symlog axis is a display transform; statistics retain the original Å values. This figure uses the mixed-backbone companion profile.

![R17: Observed-pool method complementarity ](rmsd_mixed_backbone_08_complementarity.png)

**Figure R17.** Observed-pool method complementarity using minimum RMSD as the best score. Generation counts and the available target–method cells determine the observed oracle. This figure uses the mixed-backbone companion profile.

</details>

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

### DockQ figures

All DockQ figures concern Level III RNA–protein interfaces. The conventional DockQ cutoffs are descriptive reference bands here; these analyses do not separately validate them as RNA–protein success criteria. Primary analyses retain all 144 eligible targets, with 8D97 absent for ESMFold2.

![D1: All-decoy and across-target violin distributions ](dockq_01_distributions.png)

**Figure D1.** The left panel includes all 22,000 Level III decoys. The middle and right panels show the distributions of per-target median and mean DockQ, respectively. ESMFold2 covers 143 targets; the other methods cover 144. The red line is the pooled all-decoy median (0.3208).

![D2: Mean ± SD and median ](dockq_02_mean_sd_median.png)

**Figure D2.** Mean, sample SD, and median DockQ for all decoys, target medians, and target means. Dispersion across decoys and dispersion across target summaries are distinct quantities.

![D3: Exact ECDFs ](dockq_03_ecdf.png)

**Figure D3.** Exact empirical cumulative distributions of DockQ for the three observation units. Higher scores indicate better interface agreement.

![D4: Reference-threshold coverage ](dockq_04_reference_coverage.png)

**Figure D4.** Cumulative reference coverage at DockQ ≥ 0.23, ≥ 0.49, and ≥ 0.80, for the observation units shown. Unlike the mutually exclusive categories in Table 7, these cumulative coverages overlap.

![D5: Target×method mean, SD, median and best ](dockq_05_cell_summaries.png)

**Figure D5.** Target–method mean, sample SD, median, and maximum observed DockQ across all eligible Level III targets.

![D6: Target×method all-decoy reference coverage ](dockq_06_cell_coverage.png)

**Figure D6.** Within-target, within-method fractions of decoys reaching the three cumulative DockQ reference cutoffs. Missing observations remain missing.

![D7: All-decoy target landscape ](dockq_07_target_landscape.png)

**Figure D7.** All-decoy DockQ distributions across Level III targets, separated by method and accompanied by target-level summaries.

![D8: Observed-pool complementarity ](dockq_08_complementarity.png)

**Figure D8.** Complementarity of the observed methods in reaching the best available interface score. Oracle gains depend on the actual generated decoy pools.

![D9: Across-target method correlations ](dockq_09_correlations.png)

**Figure D9.** Across-target method correlations for target medians and target means, using the 143 targets common to all methods for this comparison. Primary results retain all 144 eligible targets.

![D10: All-decoy sequence length vs DockQ ](dockq_10_length.png)

**Figure D10.** All-decoy DockQ versus reference RNA + protein sequence length for Level III complexes. Each point is one decoy; the length axis is logarithmic. All 22,000 decoys are retained without jitter or subsampling.

## Interpretation and limitations

- Decoys from the same target are related observations. Pooled distributions describe the generated pool; target-level comparisons use targets as the statistical unit.
- Method output counts differ. Best-observed scores, winner shares, and oracle comparisons reflect the actual observed generation budgets, not equal-budget comparisons.
- Level III complex TM-score and mixed-backbone metrics have a different evaluation scope from RNA-only metrics. Differences between levels also reflect differences in target composition and length.
- RMSD is calculated on matched atom pairs after a global superposition. It does not by itself penalize all unmatched residues; atom coverage should be considered when reusing these scores.
- Missingness is metric dependent. Comparisons requiring a shared population should explicitly align valid targets and decoys rather than assume that all four metrics have identical denominators.
- Score dispersion describes variation in reference similarity; it is not a direct measure of structural diversity. These analyses also do not establish QA-ranking performance or performance on unrelated future targets.


