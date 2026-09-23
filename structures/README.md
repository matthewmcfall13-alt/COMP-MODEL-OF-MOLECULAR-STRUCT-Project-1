# Protein Structures

This directory contains the experimental reference structures and the highest-ranked AlphaFold2/ColabFold predictions used in the final analysis.

## Experimental Structures

Located in:

`structures/experimental/`

| Target | Experimental PDB | Reference Chain Used |
|---|---|---|
| T1188 Easy | 8C6Z | A |
| T1106s1 Moderate | 7QIH | B |
| T1122 Hard | 8BBT | A |

## Predicted Structures

Located in:

`structures/predicted/`

The rank-1 prediction from each ColabFold run was retained.

### T1188 Easy

| MSA Condition | Rank-1 AlphaFold2 Model |
|---|---|
| Automatic MMseqs2 | model_3_seed_000 |
| Custom HMMER | model_3_seed_000 |
| Custom DeepMSA2 | model_3_seed_000 |

### T1106s1 Moderate

| MSA Condition | Rank-1 AlphaFold2 Model |
|---|---|
| Automatic MMseqs2 | model_3_seed_000 |
| Custom HMMER | model_3_seed_000 |

### T1122 Hard

| MSA Condition | Rank-1 AlphaFold2 Model |
|---|---|
| Automatic MMseqs2 | model_4_seed_000 |
| Custom HMMER | model_5_seed_000 |

## Structural Comparison

Predicted structures were mapped to the corresponding experimental sequence and compared using the same experimentally resolved residues within each target.

Mapped residues:

- T1188: **573 / 630**
- T1106s1: **71 / 122**
- T1122: **230 / 241**

For T1106s1, only the experimentally resolved region of chain B in 7QIH was available for comparison.

The structural metrics calculated from these models are available in:

`analysis/casp15_structure_comparison_metrics.csv`
