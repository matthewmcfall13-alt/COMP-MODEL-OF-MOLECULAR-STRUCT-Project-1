# Limitations

## Limited Number of Targets

This study evaluates only three CASP15 targets:

- T1188
- T1106s1
- T1122

The results should therefore be interpreted as a case study rather than a general rule describing all proteins.

## MSA Depth Is Not the Only Variable

Differences between prediction conditions cannot be attributed to MSA depth alone.

Other factors may contribute, including:

- sequence diversity
- homolog relevance
- sequence coverage
- redundancy
- alignment quality

A smaller high-quality MSA may sometimes provide more useful evolutionary information than a much larger but noisier alignment.

## Moderate Target Experimental Coverage

For T1106s1, only **71 of 122 residues** could be mapped to the experimental structure in PDB 7QIH.

Therefore, structural metrics for this target describe only the experimentally resolved and mapped region rather than the complete protein.

## Confidence Is Not Accuracy

pLDDT and pTM are AlphaFold2 confidence measures and do not directly measure agreement with an experimental structure.

For example, the Hard target demonstrates that a prediction can have slightly higher confidence without having higher experimental structural accuracy.

## HMMER Filtering

The T1188 JackHMMER search initially produced a very large number of sequences, including many short fragments.

The final custom HMMER MSA was filtered using a minimum 70% target-coverage criterion.

This filtering step changes both the size and composition of the alignment and should be considered when comparing it with the automatic MMseqs2 and DeepMSA2 alignments.

## Missing Original Notebooks

The original prediction notebooks were not retained for:

- T1106s1 Automatic MMseqs2
- T1122 Automatic MMseqs2

However, their MSA inputs, predicted rank-1 structures, and final analysis results are retained in the repository.

## Structural Metrics

Different structural metrics emphasize different properties.

For example:

- TM-score emphasizes overall fold similarity.
- lDDT emphasizes preservation of local structural environments.
- RMSD is sensitive to displacement and protein size.

Consequently, one prediction may perform better according to one metric but worse according to another.

## Overall Interpretation

The results support the idea that the usefulness of an MSA depends not only on the number of sequences but also on the quality and relevance of the evolutionary information it contains.

However, a larger benchmark containing many more proteins would be required to test this relationship statistically.
