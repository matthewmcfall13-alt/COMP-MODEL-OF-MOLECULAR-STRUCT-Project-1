# Figure captions and result notes for Allen

**Alignment depth and diversity (`alignment_depth_and_diversity.png`).** All supplied records versus distinct query-coordinate sequence strings for each MSA condition. Counts include the query; exact duplicates retain their multiplicity in the full-row statistics. The logarithmic axis displays the large depth range without treating sequence count as a direct accuracy measure.

**Alignment profiles (`alignments/<target>/<condition>/alignment_profiles.png`).** Per-position non-gap coverage and gap fraction use all records. Canonical-residue entropy and consensus exclude gaps/ambiguous residues; non-query profiles exclude row zero. Identity/coverage distributions describe all eligible rows after the query. Query-only or duplicate-query conditions contain no distinct homolog-sequence information.

**Alignment overview (`alignments/<target>/<condition>/alignment_overview.png`).** Query-coordinate view of up to 80 deterministically spaced file-order rows. Colors distinguish gaps, noncanonical residues, mismatches and query identity. This is a display subset; every input row contributes to numerical summaries.

**Structure profiles (`structures/<target>/per_residue_error_confidence_coverage.png`).** Cα deviations after independent rigid fits on one shared experimental mask, predicted per-residue pLDDT, and coordinate coverage. Missing experimental positions are unscored. pLDDT describes model confidence rather than measured structural accuracy.

## Target-specific notes

- **T1188, Custom HMMER versus automatic:** reported changes are 0.000 TM-score, -0.60 Cα-lDDT percentage points, and 0.01 Å RMSD (custom minus automatic). Higher TM/lDDT and lower RMSD indicate better agreement on the stated mask.
- **T1188, Custom DeepMSA2 versus automatic:** reported changes are 0.001 TM-score, -0.25 Cα-lDDT percentage points, and -0.03 Å RMSD (custom minus automatic). Higher TM/lDDT and lower RMSD indicate better agreement on the stated mask.
- **T1106s1, Custom HMMER versus automatic:** reported changes are -0.029 TM-score, 2.20 Cα-lDDT percentage points, and 2.96 Å RMSD (custom minus automatic). Higher TM/lDDT and lower RMSD indicate better agreement on the stated mask.
- **T1122, Custom HMMER versus automatic:** reported changes are 0.070 TM-score, 5.46 Cα-lDDT percentage points, and -0.87 Å RMSD (custom minus automatic). Higher TM/lDDT and lower RMSD indicate better agreement on the stated mask.
- **T1122 has no distinct sequence beyond its query in the supplied alignments.** Three automatic records are duplicate query sequences; the HMMER alignment is query-only. Their counts do not provide three independent homologs.
- **T1106s1 coverage:** 71/122 target residues are scored. Its structural results describe the experimentally mapped segment, not the entire target.

Nancy's existing overview/accuracy figures remain in `../../figures/`; cite `../../analysis/casp15_structure_comparison_metrics.csv` for their reported values. The new figures supplement those results with alignment composition and residue-level evidence. Main comparisons are HMMER versus automatic; the T1188 DeepMSA2 result is an additional comparison.
