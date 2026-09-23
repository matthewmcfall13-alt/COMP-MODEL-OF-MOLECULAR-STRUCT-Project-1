# Local analysis of the completed AF2/ColabFold study

The current study uses T1188, T1106s1 and T1122. Each has automatic MMseqs2 and custom JackHMMER inputs and predictions; T1188 additionally has a DeepMSA2 condition. This accepted AF2/ColabFold design supersedes the earlier AF3 submission proposal. Nancy's [methods](methods.md), [limitations](limitations.md) and [analysis summary](../analysis/analysis_summary.md) document the experiments and reported results.

Predictions used ColabFold v1.6.3. The automatic mode was `mmseqs2_uniref_env`; custom JackHMMER used one iteration against UniProt 2025_01. For T1188, custom hits were filtered to at least 70% target coverage. The analysis retains these as the team's experimental methods; additional prediction jobs or independent verification of every reported value are not prerequisites for completing the project.

## Reproduction and output locations

Follow [MATTHEW_WORKFLOW.md](../MATTHEW_WORKFLOW.md) for Python 3.10+ setup. Run all cells of the code-only `Project1_Server_Workflow.ipynb`, or run:

```powershell
.\.venv\Scripts\python.exe -m project1 --mode uploaded
```

The `uploaded_analysis` configuration maps the source files explicitly. Local analysis reads Nancy's inputs and writes its own outputs. Detailed snapshots are retained under `outputs/analyses/uploaded/`; stable tables and figures are published under [results/matthew_analysis/](../results/matthew_analysis/REPORT.md). The [results guide](../RESULTS_GUIDE.md) identifies figures for Allen and Yu Wu.

## Alignment calculations

Original A3M bytes and row order are preserved. The separate query-coordinate analysis removes lowercase insertion letters and insertion-gap dots. ColabFold's single-chain length header is recognized and checked. The query must match the target sequence and every analyzed row must have the target's match-column width.

- **Raw depth** includes the query. Exact duplicate counts compare complete query-coordinate strings, including gap positions. Distinct non-query sequences are unique strings different from the query; this is not a count of independent evolutionary observations.
- **Coverage, gaps and identity** describe alignment occupancy and similarity. Canonical amino-acid identities exclude gaps and ambiguous symbols from the eligible denominator.
- **BLOSUM62 sum-of-pairs** is normalized per eligible unordered canonical residue pair. It has matrix-score units, not a universal probability or alignment-quality scale.
- **Conservation** uses canonical-residue Shannon entropy, excluding gaps and ambiguous symbols. Normalized conservation is `1 - entropy/log2(20)`. Positions without eligible residues have missing statistics.
- **Neff**, where available, is an exact bounded calculation with the definition in [METRICS.md](METRICS.md). It is not computed automatically above 500 rows. Nancy's DeepMSA-reported Nf is retained separately; it is not relabeled as locally computed Neff.

All-row statistics use the complete alignments. Alignment overview figures use a deterministic row subset only for display and label the displayed and total counts. T1122 automatic contains three identical query records, and its custom HMMER contains the query alone. Conservation in these files does not demonstrate homolog conservation because neither has non-query sequence diversity.

## Coordinates, confidence and reported scores

Reference chains are 8C6Z A for T1188, 7QIH B for T1106s1 and 8BBT A for T1122. Sequence-based mapping gives a shared comparison mask within each target: **573/630**, **71/122** and **230/241** residues respectively. Residue mapping tables retain the scored and missing positions. These scores describe the mapped experimental regions, especially the partial T1106s1 chain.

Local Cα RMSD uses a Kabsch rigid fit without reflection or distance pruning. Common-mask Cα-lDDT uses reference contacts below 15 Å and distance-error thresholds of 0.5, 1, 2 and 4 Å, aggregated over eligible contacts. It is not all-atom lDDT. The local mean Cα pLDDT comes from predicted ColabFold PDB B-factor fields across the model; experimental B-factors are never interpreted as pLDDT.

The combined table retains Nancy's TM-score, pTM and other reported metrics with their source and distinguishes them from locally calculated values. The team's TM-scores use mapped-residue normalization as documented in its analysis summary. No one-fit substitute is presented as an optimized TM-score. RMSD100 is a separate size-normalized quantity and must not be labeled raw RMSD. The team's reported values remain usable without an additional external scoring installation.

## Interpretation

Paired differences are custom minus automatic: negative RMSD differences and positive lDDT/TM-score differences favor custom. pLDDT and pTM are confidence measures rather than experimental accuracy. The primary comparison is HMMER versus automatic for all three targets; T1188 DeepMSA2 is an additional comparison.

MSA depth, diversity, filtering and target biology vary together. The study supports target-specific observations, not a general causal relationship or statistical superiority across proteins. Synthetic tests validate software behavior separately from the real study results. The local workflow runs no prediction models or sequence-search services and does not alter original team files.
