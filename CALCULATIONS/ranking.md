# Case ranks and aggregate statistics

## Generic population benchmark

1. For each of the 16 spectrum-by-budget keys, obtain the full-window \(P_1\) relative \(L_2\) error. Use `Pz_relL2_vs_HEOM` in the [population summary](../DATA/mDA/generic_zero_T/population_summary.csv) for Linear, Log, MDM-c, BSDO, LMD-A, LMD-B and LMD-S. Use the corresponding field in the [canonical ID summary](../DATA/mDA/generic_zero_T/canonical_id_generic_summary.csv) for ID. The trajectory CSVs provide independent checks of those values.
2. Rank the eight errors in ascending order inside each key. The archived canonical ranking used a sorted ordinal rank with method order Linear, Log, MDM, BSDO, LMD-A, LMD-B, LMD-S, ID for exact ties. Rank 1 is the smallest error.
3. For each method, calculate the arithmetic mean and median of the case ranks, the number of rank-1 cases, the maximum error, and \(\exp(\operatorname{mean}(\log e))\) over the positive case errors. The original geometric-mean implementation floors errors at the smallest positive normalized double before taking logs. Sort the aggregate by mean rank and then geometric-mean error.
4. Compare the output with [the archived case ranks](../DATA/mDA/generic_zero_T/generic_method_ranking_by_case_with_canonical_id.csv) and [aggregate](../DATA/mDA/generic_zero_T/generic_method_ranking_with_canonical_id.csv).

## BCF–population Spearman correlation

Inside each spectrum–budget key, join the full BCF scores from the two generic BCF batches with the canonical ID BCF score. In the same method order, pair them with the HEOM population errors. Assign average ranks to ties separately for the BCF and population lists, then take the Pearson correlation of those two rank vectors. Repeat after removing ID before ranking. The summary gives the arithmetic mean, minimum, maximum and number strictly below zero for each list of 16 coefficients.

## BChl distances

The representative fixed-trajectory table ranks eight methods by the full relative distance to the fixed BSDO–mD₂ curve within each budget, then aggregates across 10, 20, 30 and 40 modes. Use the canonical ID distance to that same trajectory. The separate 18-case grid ranks only LMD-A/B/S within each linewidth–coupling–budget key against its matched fixed trajectory. In a sorted list, the sensitivity script groups a value \(x\) with the smallest value \(x_0\) in the current tie group if \(|x-x_0|\le\max(10^{-14},10^{-12}\max(|x|,|x_0|))\). It assigns the group the mean of its ordinal ranks. For the linewidth–coupling summary, compute the geometric mean over the two included budgets per method. Preserve the two different grouping keys when making comparisons.

## November depth-margin screen

Each method's RMS error \(e_i\) is paired with the corresponding RMS change \(d_i\) between adjacent HEOM depths. Compute \(L_i=\max(0,e_i-d_i)\) and \(U_i=e_i+d_i\). The saved screen assigns \(1+\#\{j:U_j<L_i\}\) within a complete case and time window. Methods can share that rank. This calculation is a screen against the recorded adjacent-depth margin.
