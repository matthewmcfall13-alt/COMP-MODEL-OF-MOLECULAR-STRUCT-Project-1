# Matthew's alignment and structure results

**Status: complete for the seven supplied study conditions.**

Source Git revision: `9ad1fb66039963b66a82680c3bf70a9f6c9e7e36`. Generated: 2026-09-23T13:46:59.963910+00:00. 7/7 condition rows combine an analyzed alignment, mapped coordinates and the team's result table.

The accepted study uses Nancy's T1188, T1106s1 and T1122 targets with AlphaFold2/ColabFold. The main comparison is automatic MMseqs2 versus custom HMMER; T1188 additionally compares DeepMSA2. Nancy's structural results are accepted team results and retained in the `team_` columns. Matthew's alignment statistics and coordinate recalculations appear separately; neither replaces the original measurements.

[Team methods](../../docs/methods.md) · [Team limitations](../../docs/limitations.md) · [Original result table](../../analysis/casp15_structure_comparison_metrics.csv) · [Original structure figures](../../figures/)

## Combined results

| Target | Condition | Rows / distinct sequences | TM-score | Cα-lDDT (%) | RMSD (Å) | Available inputs |
|---|---|---:|---:|---:|---:|---|
| T1188 | Automatic MMseqs2 | 11337 / 11295 | 0.990 | 96.88 | 0.87 | COMPLETE |
| T1188 | Custom HMMER | 213 / 188 | 0.990 | 96.28 | 0.88 | COMPLETE |
| T1188 | Custom DeepMSA2 | 25136 / 25136 | 0.991 | 96.63 | 0.84 | COMPLETE |
| T1106s1 | Automatic MMseqs2 | 113 / 105 | 0.674 | 86.28 | 3.82 | COMPLETE |
| T1106s1 | Custom HMMER | 130 / 100 | 0.645 | 88.48 | 6.78 | COMPLETE |
| T1122 | Automatic MMseqs2 | 3 / 1 | 0.426 | 50.44 | 20.82 | COMPLETE |
| T1122 | Custom HMMER | 1 / 1 | 0.496 | 55.90 | 19.95 | COMPLETE |

Structural values in this table are Nancy's accepted reported values; complete precision for recalculated values is in `combined_results.csv`. Blank values remain missing, never zero.

## What the comparison shows

- **T1188, Custom HMMER versus automatic:** reported changes are 0.000 TM-score, -0.60 Cα-lDDT percentage points, and 0.01 Å RMSD (custom minus automatic). Higher TM/lDDT and lower RMSD indicate better agreement on the stated mask.
- **T1188, Custom DeepMSA2 versus automatic:** reported changes are 0.001 TM-score, -0.25 Cα-lDDT percentage points, and -0.03 Å RMSD (custom minus automatic). Higher TM/lDDT and lower RMSD indicate better agreement on the stated mask.
- **T1106s1, Custom HMMER versus automatic:** reported changes are -0.029 TM-score, 2.20 Cα-lDDT percentage points, and 2.96 Å RMSD (custom minus automatic). Higher TM/lDDT and lower RMSD indicate better agreement on the stated mask.
- **T1122, Custom HMMER versus automatic:** reported changes are 0.070 TM-score, 5.46 Cα-lDDT percentage points, and -0.87 Å RMSD (custom minus automatic). Higher TM/lDDT and lower RMSD indicate better agreement on the stated mask.
- **T1122 has no distinct sequence beyond its query in the supplied alignments.** Three automatic records are duplicate query sequences; the HMMER alignment is query-only. Their counts do not provide three independent homologs.
- **T1106s1 coverage:** 71/122 target residues are scored. Its structural results describe the experimentally mapped segment, not the entire target.

These are descriptive within-target comparisons. Search method, filtering, redundancy, coverage and model selection vary alongside depth; this three-target case study does not isolate a causal effect of sequence count. Small numerical differences are not evidence of statistical significance.

## Files and metric conventions

- `combined_results.csv`: joined alignment statistics, accepted team structural/confidence metrics and separate coordinate recalculations for every planned condition.
- `paired_custom_minus_auto.csv`: signed differences. Negative RMSD/RMSD100 and positive TM-score/lDDT favor custom; pLDDT/pTM measure confidence, not experimental accuracy.
- `source_inventory.csv`: repository-relative source paths and original SHA-256 hashes. Team inputs are read only.
- `alignment_depth_and_diversity.png` and `alignments/`: full-row coverage, gaps, query identity, redundancy, BLOSUM62 sum-of-pairs, entropy, and labeled alignment overviews. Overview rows are a deterministic display subset; numerical statistics use all rows.
- `structures/`: shared residue masks, residue mappings, per-residue error/confidence/coverage figures and recalculated comparisons.
- `FIGURE_CAPTIONS.md`: text Allen can use with these figures independently of a slide deck.

Counts include the query. Distinct sequences are exact insertion-stripped query-coordinate strings, including gap positions; they are not a count of independent evolutionary observations. Entropy excludes gaps/noncanonical symbols. BLOSUM62 sum-of-pairs excludes pairs containing gaps/noncanonical residues and is normalized by eligible pair count. Exact local Neff is only calculated within the documented 500-row limit; a blank value for a large alignment is unavailable, not zero, and is distinct from DeepMSA's reported Nf.

Coordinate recalculations use the same mapped Cα mask for all conditions of each target. Kabsch RMSD uses no outlier pruning. Cα-lDDT uses reference contacts under 15 Å and strict distance-error thresholds 0.5/1/2/4 Å. Whole-model mean pLDDT uses one value per mapped model residue. Nancy's TM-score normalization is the mapped experimental count; confidence and accuracy remain separate.

## Reproduction

From the repository, run `python -m project1 --mode uploaded` or all code cells in `Project1_Server_Workflow.ipynb`. Publication is enabled by `uploaded_analysis.publish_results` in `config.json`. Only files owned by this generated folder are refreshed; original team folders are preserved.
