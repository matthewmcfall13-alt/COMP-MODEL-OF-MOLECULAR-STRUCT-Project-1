# Methods

## Study Design

Three CASP15 single-chain targets were selected to represent different levels of homolog availability:

- **Easy:** T1188 / PDB 8C6Z
- **Moderate:** T1106s1 / PDB 7QIH
- **Hard:** T1122 / PDB 8BBT

The main comparison was:

**Automatic ColabFold MMseqs2 MSA vs Custom JackHMMER MSA**

For the Easy target, a third condition using **DeepMSA2** was also evaluated.

## MSA Construction

### Automatic MSA

Automatic alignments were generated through ColabFold using:

`msa_mode = mmseqs2_uniref_env`

### Custom HMMER MSA

Custom alignments were generated using JackHMMER against:

**UniProt 2025_01**

Only **Iteration 1** was used for the final analysis to avoid profile drift during later iterations.

Final custom HMMER MSA depths were:

| Target | MSA Depth |
|---|---:|
| T1188 Easy | 213 |
| T1106s1 Moderate | 130 |
| T1122 Hard | 1 |

For T1188, the raw HMMER search contained many short fragments. Hits were therefore filtered to retain sequences covering at least 70% of the target sequence.

### DeepMSA2

DeepMSA2 was additionally used for T1188.

The resulting alignment contained:

**25,136 sequences**

with an effective depth:

**Nf ≈ 285.7**

## Structure Prediction

Protein structures were predicted using:

**AlphaFold2 through ColabFold v1.6.3**

Each MSA condition was submitted independently.

The highest-ranked model from each prediction was used for structural comparison.

## Experimental References

The experimental structures were:

| Target | PDB | Method | Resolution |
|---|---|---|---:|
| T1188 | 8C6Z | X-ray crystallography | 1.85 Å |
| T1106s1 | 7QIH | X-ray crystallography | 1.92 Å |
| T1122 | 8BBT | X-ray crystallography | 1.69 Å |

## Evaluation Metrics

Predicted structures were compared against the corresponding experimental structure using:

- pLDDT
- pTM
- TM-score
- Cα-lDDT
- Cα RMSD
- RMSD100

pLDDT and pTM were treated as model-confidence measures.

TM-score, Cα-lDDT, RMSD, and RMSD100 were used as experimental-structure accuracy measures.

## Structural Mapping

Only residues that could be mapped between the prediction and experimental structure were used for structural scoring.

Mapped residues:

| Target | Mapped Experimental Residues |
|---|---:|
| T1188 | 573 / 630 |
| T1106s1 | 71 / 122 |
| T1122 | 230 / 241 |

For T1106s1, the comparison used the resolved portion of **chain B of PDB 7QIH**.

Because only 71 of 122 target residues were experimentally resolved and mapped, results for this target should be interpreted with additional caution.

## Structural Superposition

Cα atoms were superimposed using least-squares Kabsch alignment.

Cα-lDDT was calculated using a 15 Å reference-neighbor cutoff and distance tolerances of:

- 0.5 Å
- 1 Å
- 2 Å
- 4 Å

RMSD100 was used to reduce the dependence of raw RMSD on protein length.

## Interpretation

The experiment was designed as a three-target case study.

Therefore, the results illustrate how MSA depth and homolog quality may influence AlphaFold2 predictions across these selected targets, but they should not be interpreted as a universal relationship for all proteins.
