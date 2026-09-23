# Target Selection

## Objective

Three CASP15 single-chain targets were selected to represent different levels of evolutionary information available for multiple sequence alignment (MSA) construction:

- **Easy:** deep homolog availability
- **Moderate:** intermediate homolog availability
- **Hard:** very shallow or essentially absent homolog availability

A second selection criterion was the availability of a high-quality experimental structure determined by X-ray crystallography at a resolution of ≤ 2.0 Å.

## Selected Targets

| Difficulty | CASP15 Target | PDB | Length | X-ray Resolution | Homolog Availability |
|---|---|---|---:|---:|---|
| Easy | T1188 | 8C6Z | 630 aa | 1.85 Å | Deep |
| Moderate | T1106s1 | 7QIH | 122 aa | 1.92 Å | Intermediate |
| Hard | T1122 | 8BBT | 241 aa | 1.69 Å | Very shallow |

## Easy Target — T1188 / 8C6Z

T1188 was selected as the Easy target because sequence searches identified a large homologous protein family.

DeepMSA2 produced an MSA containing **25,136 sequences** with an effective depth (Nf) of approximately **285.7**.

A one-iteration JackHMMER search also identified a very large number of homologous sequences. After filtering the alignment to retain sequences covering at least 70% of the 630-residue target, **212 homologs** were retained for the custom HMMER MSA.

The experimental reference structure, **PDB 8C6Z**, was determined by X-ray crystallography at **1.85 Å resolution**.

These observations support classification of T1188 as a protein with a **deep and information-rich sequence family**.

## Moderate Target — T1106s1 / 7QIH

T1106s1 was selected as the Moderate target because its homolog availability was substantially lower than T1188 but clearly greater than T1122.

A one-iteration JackHMMER search produced an alignment containing **130 sequences**. Most identified sequences were annotated as YscX, AscX, or related type III secretion proteins.

The experimental reference structure, **PDB 7QIH**, was determined by X-ray crystallography at **1.92 Å resolution**.

T1106s1 therefore represents an **intermediate-depth protein family** between the Easy and Hard targets.

## Hard Target — T1122 / 8BBT

T1122 was selected as the Hard target because multiple sequence-search approaches found essentially no additional homologous sequence information.

DeepMSA2 returned only **one sequence**, corresponding to the query itself.

A one-iteration JackHMMER search likewise returned only the target sequence, with no additional usable homologs.

The experimental reference structure, **PDB 8BBT**, was determined by X-ray crystallography at **1.69 Å resolution**.

The agreement between independent search methods supports classification of T1122 as an **extremely shallow protein family with essentially no usable evolutionary information**.

## Selection Rationale

The final target set therefore spans a broad range of MSA depth while maintaining high-quality experimental reference structures:

**T1188 (Easy/deep) → T1106s1 (Moderate/intermediate) → T1122 (Hard/shallow)**

This design allows the effect of homolog availability and MSA construction strategy on AlphaFold2 structure prediction to be investigated across three distinct evolutionary-information regimes.
