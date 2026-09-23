# T1188 Easy — MSA Inputs

Target: **T1188**  
Experimental structure: **PDB 8C6Z**  
Difficulty class: **Easy / deep homolog availability**  
Target length: **630 aa**

This folder contains the three MSA conditions evaluated for T1188.

## Automatic MMseqs2

File:

`T1188_Easy_Automatic_MMseqs2.a3m`

Generated automatically by ColabFold using the `mmseqs2_uniref_env` search mode.

MSA depth used for prediction:

**11,337 sequences**

## Custom HMMER

File:

`T1188_Easy_Custom_HMMER.a3m`

Generated using one JackHMMER iteration against UniProt (2025_01).

The raw HMMER result contained 59,835 sequences and was strongly enriched in chitinase-family homologs.

Because many hits were short fragments, the alignment was filtered to retain sequences covering at least 70% of the 630-residue target.

Final MSA depth:

**213 sequences total**
- 1 query sequence
- 212 retained homologs

## Custom DeepMSA2

File:

`T1188_Easy_Custom_DeepMSA2.a3m`

Generated using DeepMSA2.

Final MSA depth:

**25,136 sequences**

DeepMSA2 reported an effective alignment depth:

**Nf ≈ 285.7**

## Purpose

The three alignments allow T1188 to be used as a deep-family test case for comparing:

- automatic MSA generation
- a filtered HMMER MSA
- a much deeper DeepMSA2 MSA

Despite the large difference in raw MSA depth, all three conditions produced very similar high-accuracy AlphaFold2 structures.
