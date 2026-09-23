# T1122 Hard — MSA Inputs

Target: **T1122**  
Experimental structure: **PDB 8BBT**  
Difficulty class: **Hard / extremely shallow homolog availability**  
Target length: **241 aa**

This folder contains the two MSA conditions evaluated for T1122.

## Automatic MMseqs2

File:

`T1122_Hard_Automatic_MMseqs2.a3m`

Generated automatically by ColabFold using the `mmseqs2_uniref_env` search mode.

MSA depth used for prediction:

**3 sequences**

## Custom HMMER

File:

`T1122_Hard_Custom_HMMER.a3m`

Generated from a one-iteration JackHMMER search against UniProt (2025_01).

JackHMMER returned only the target sequence, with no additional usable homologs.

Final MSA depth:

**1 sequence**

This is therefore effectively a single-sequence prediction condition.

## Independent Evidence for Shallow Homology

DeepMSA2 also returned only **1 sequence**, corresponding to the query itself.

The agreement between JackHMMER and DeepMSA2 supports the classification of T1122 as an extremely shallow target with essentially no usable evolutionary information.

## Purpose

T1122 serves as the low-information endpoint of the Easy–Moderate–Hard target series.

Its Auto and Custom predictions test AlphaFold2 performance when little or no evolutionary information is available.
