# Figures

This directory contains the quantitative comparison plots and structural overlays used to summarize the CASP15 MSA experiment.

## Quantitative Figures

### `figure_tm_score.png`

Comparison of TM-score across MSA conditions.

TM-score measures overall fold similarity between the predicted and experimental structures. Values closer to 1 indicate greater global structural agreement.

### `figure_lddt_ca.png`

Comparison of Cα-lDDT across MSA conditions.

Cα-lDDT measures preservation of local structural environments without requiring a global superposition.

### `figure_rmsd100.png`

Comparison of RMSD100 across prediction conditions.

RMSD100 is a size-normalized version of RMSD that facilitates comparison between proteins of different lengths. Lower values indicate better agreement with the experimental structure.

### `figure_msa_depth.png`

Comparison of the number of sequences used in each MSA condition.

The large differences in MSA depth illustrate that alignment size alone does not determine prediction accuracy.

### `figure_easy_three_way_tm_score.png`

Three-way TM-score comparison for the Easy target T1188:

- Automatic MMseqs2
- Custom HMMER
- Custom DeepMSA2

Despite large differences in MSA depth, all three approaches produced very similar high-accuracy structures.

## Structural Overlays

### `overlay_T1188_easy.png`

Cα structural comparison for the Easy target T1188.

### `overlay_T1106s1_moderate.png`

Cα structural comparison for the Moderate target T1106s1.

### `overlay_T1122_hard.png`

Cα structural comparison for the Hard target T1122.

These overlays were generated from the corresponding experimental and predicted PDB coordinates after structural superposition.

They are intended as simplified backbone visualizations rather than publication-style molecular ribbon renderings.

## Numerical Data

The values underlying these figures are available in:

`analysis/casp15_structure_comparison_metrics.csv`
