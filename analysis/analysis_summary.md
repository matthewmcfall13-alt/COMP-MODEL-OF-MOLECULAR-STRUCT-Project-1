# CASP15 MSA-depth / AlphaFold2–ColabFold analysis

## Research question
How does MSA depth and MSA-generation strategy affect AlphaFold2 structure prediction across proteins with deep, intermediate, and shallow homolog availability?

## Experimental targets
| Difficulty | CASP target | Experimental PDB | X-ray resolution | Residues used for experimental scoring |
|---|---|---|---:|---|
| Easy | T1188 | 8C6Z | 1.85 Å | 573 / 630 (target positions 25–597) |
| Moderate | T1106s1 | 7QIH chain B | 1.92 Å | 71 / 122 (target positions 50–120) |
| Hard | T1122 | 8BBT chain A | 1.69 Å | 230 / 241 (target positions 4–237; four unresolved internal residues plus terminal residues excluded) |

All Auto/Custom comparisons for a given target use the exact same experimentally resolved residue mapping.

## Main Auto vs Custom-HMMER results
| Difficulty | Condition | MSA depth | pLDDT | pTM | TM-score | Cα-lDDT | Cα RMSD (Å) | RMSD100 (Å) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Easy | Automatic MMseqs2 | 11,337 | 90.66 | 0.885 | 0.990 | 96.88 | 0.87 | 0.46 |
| Easy | Custom HMMER | 213 | 89.72 | 0.878 | 0.990 | 96.28 | 0.88 | 0.47 |
| Moderate | Automatic MMseqs2 | 113 | 70.69 | 0.413 | 0.674 | 86.28 | 3.82 | 4.60 |
| Moderate | Custom HMMER | 130 | 72.63 | 0.422 | 0.645 | 88.48 | 6.78 | 8.18 |
| Hard | Automatic MMseqs2 | 3 | 45.42 | 0.493 | 0.426 | 50.44 | 20.82 | 14.70 |
| Hard | Custom HMMER | 1 | 45.82 | 0.464 | 0.496 | 55.91 | 19.95 | 14.08 |

## Easy-target DeepMSA2 bonus comparison
| T1188 condition | MSA depth | pLDDT | pTM | TM-score | Cα-lDDT | Cα RMSD (Å) | RMSD100 (Å) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Automatic MMseqs2 | 11,337 | 90.66 | 0.885 | 0.990 | 96.88 | 0.87 | 0.46 |
| Custom HMMER | 213 | 89.72 | 0.878 | 0.990 | 96.28 | 0.88 | 0.47 |
| Custom DeepMSA2 | 25,136 | 90.86 | 0.885 | 0.991 | 96.63 | 0.84 | 0.45 |

## Interpretation
1. **Easy / deep family:** all three MSA strategies produced essentially the same high-quality structure. Increasing the MSA from 213 filtered HMMER sequences to 25,136 DeepMSA2 sequences produced only a very small accuracy change, consistent with a saturation / diminishing-returns regime for this target.
2. **Moderate / intermediate family:** the result is metric-dependent. Custom HMMER has higher pLDDT, pTM and Cα-lDDT, indicating somewhat better confidence/local geometry, while Automatic MMseqs2 has better TM-score and much lower global RMSD/RMSD100. This suggests the custom model preserves local neighborhoods well but differs more in the global placement of the experimentally resolved segment.
3. **Hard / shallow family:** both predictions are low-confidence and globally poor. The query-only HMMER condition nevertheless gives a higher experimental TM-score and lDDT and slightly lower RMSD than the 3-sequence automatic MSA. The difference should be treated as target-specific because neither MSA contains meaningful evolutionary depth.
4. **MSA depth alone is not sufficient:** 213 carefully filtered Easy-target HMMER sequences perform almost identically to 11,337 Automatic and 25,136 DeepMSA2 sequences. Homolog relevance/quality and the information content of the family matter in addition to raw sequence count.
5. **Confidence is not experimental accuracy:** for the Hard target, Automatic has higher pTM than Custom HMMER, but Custom HMMER agrees better with 8BBT by TM-score, lDDT and RMSD.

## Structural scoring method
- Rank-1 ColabFold model was used for every condition.
- Predictions and experimental structures were mapped through the known target sequence.
- Only residues resolved in the experimental PDB were scored, and the same residue mask was used for all conditions within each target.
- Cα RMSD: least-squares Kabsch superposition on all mapped Cα atoms.
- Cα-lDDT: reference-neighbor cutoff 15 Å; distance tolerances 0.5, 1, 2 and 4 Å.
- TM-score: standard TM-score distance function, normalized to the number of mapped experimental residues, with the fixed target-to-PDB residue correspondence and TM-score-optimized superposition.
- RMSD100: Carugo–Pongor size normalization: RMSD / [1 + ln(sqrt(N/100))].

## Important limitation
Only three targets were studied, so this is a controlled case study rather than evidence of a universal relationship between MSA depth and AlphaFold2 accuracy. T1106s1 also has only 71 experimentally resolved target residues in 7QIH chain B, so its global metrics should be interpreted specifically for that resolved segment.
