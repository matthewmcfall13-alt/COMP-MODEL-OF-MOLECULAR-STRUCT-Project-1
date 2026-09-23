# T1106s1 Moderate — MSA Inputs

Target: **T1106s1**  
Experimental structure: **PDB 7QIH**  
Difficulty class: **Moderate / intermediate homolog availability**  
Target length: **122 aa**

This folder contains the two MSA conditions evaluated for T1106s1.

## Automatic MMseqs2

File:

`T1106s1_Moderate_Automatic_MMseqs2.a3m`

Generated automatically by ColabFold using the `mmseqs2_uniref_env` search mode.

MSA depth used for prediction:

**113 sequences**

## Custom HMMER

File:

`T1106s1_Moderate_Custom_HMMER.a3m`

Generated using one JackHMMER iteration against UniProt (2025_01).

The downloaded Iteration-1 alignment contained **130 sequences**.

Most homologs were annotated as:

- YscX family type III secretion proteins
- AscX proteins
- related type III secretion-system proteins

The exact 122-residue T1106s1 query was placed first in the final A3M.

## Purpose

T1106s1 represents an intermediate-depth sequence family, providing a comparison point between the deep T1188 family and the extremely shallow T1122 family.

The Automatic and Custom HMMER alignments were used independently as inputs to AlphaFold2 through ColabFold.
