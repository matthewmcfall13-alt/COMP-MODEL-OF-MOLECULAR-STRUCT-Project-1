# Implemented metric definitions and limits

The final accepted AF2 study uses the alignment formulas below and the supplied-model analysis described in [UPLOADED_ANALYSIS.md](UPLOADED_ANALYSIS.md). Nancy's reported TM-score/pTM values remain the team results. Sections about server-request validation, AUTO/DEEPMSA imports, external scoring wrappers and submission preparation describe the retained historical AF3 workflow, not additional requirements for the completed AF2 project.

`project1/metrics.py` performs CPU-only local analysis. A generated request, parsed
fixture, or passing local test is not a server result. The initial FASTA-only run
has no measured structural accuracy. `metric_summary.json` lists availability;
failed or unavailable scores remain JSON `null`, never zero.

## Alignment coordinates and denominators

The original A3M text is retained and submitted unchanged after validation. A
separate query-coordinate analysis copy removes lowercase insertion residues and
`.` insertion-gap characters. Uppercase residues and `-` match-column gaps remain
in their original order. The first sequence must exactly equal the submitted
query, without gaps or lowercase insertions. Every analysis row has query length
L. No alignment is truncated and the project targets are never aligned together.

All counts include the query unless the field explicitly says otherwise:

- Sequence count N includes the query; homolog-row count is N minus one.
- Per-position non-gap coverage is non-`-` rows divided by N. Gap fraction is gap
  rows divided by N. Canonical coverage counts only the 20 standard amino acids.
- Identity to query for each homolog is identical canonical residues divided by
  that homolog's canonical-residue count in query columns. Gaps and ambiguous
  residues do not enter the denominator. An all-gap row has missing identity.
  The reported mean excludes the query and rows with missing identity.
- Duplicates are exact equal query-coordinate strings, including gaps; insertion
  differences are not considered. Duplicate fraction is (N minus unique rows)/N.

The BLOSUM62 sum-of-pairs uses the Biopython BLOSUM62 matrix. For column j with
canonical residue counts n(j,a), its score is

`S_j = sum_a choose(n(j,a),2)*B(a,a) + sum_(a<b) n(j,a)*n(j,b)*B(a,b)`.

The eligible-pair denominator is
`D = sum_j choose(sum_a n(j,a),2)`. The normalized score is `sum_j S_j / D`.
The query participates. All gap-containing and noncanonical-containing pairs are
excluded, with no gap penalty. When D is zero, the normalized score is missing.
This is a matrix score per eligible residue pair, not a probability or a 0–1
quality scale. Column counts avoid a quadratic loop over all sequence pairs.
The course sequence-alignment lecture pp.102–103 defines the unnormalized pair
sum; normalization and exclusion rules here are explicit analysis choices.

Exact Neff is optional and bounded to at most 500 rows by default. Two rows are
neighbors when (1) identity on their jointly canonical columns is at least 0.8,
and (2) jointly canonical column count divided by L is at least 0.5. The query is
included. Each row contributes `1 / neighbor_count`; a row with no eligible
neighbors contributes zero. Self counts only when the coverage rule passes.
Neff is the sum of these contributions. Above the limit, Neff is unavailable;
there is no automatic subsampling or approximation. These thresholds can be
changed in the Python function, and the actual definition is recorded.

No search E-values or p-values are inferred from A3M text. User-supplied search
metadata is retained separately. Depth, redundancy and SP score describe the
alignment; none alone demonstrates biological correctness or prediction quality.

## Reference assignment and residue mapping

Scoring requires a separately confirmed experimental reference with recorded
assignment evidence. Sequence similarity is not confirmation of a CASP mapping.
Missing references do not prevent submission preparation. The workflow requires
both imports to pass target/control verification and uses server-selected primary
models before consulting the reference.

Biopython reads PDB/mmCIF coordinates using author chain IDs, residue numbers and
insertion codes. Protein Cα atoms are retained; residues lacking Cα are recorded
as missing. The selected alternate location and occupancy remain available in
the parsed record. MSE maps to methionine for sequence matching. Multiple
coordinate models require an explicit model selection in the metrics API.

Exact target/coordinate sequences map directly. Otherwise a global sequence
alignment uses match +2, mismatch −3, gap-open −5, gap-extension −0.5. Multiple
optimal alignments or mismatched mapped residues block automatic mapping.
Explicit mappings can resolve a known construct/missing-region ambiguity:
`{ "1": "42", "2": "42A" }` maps target positions to author residue IDs in the
declared chain. Explicit mappings still require unique, sequence-identical pairs.
Omitted target positions remain unscored. Unmapped coordinate residues and missing
Cα records are retained in provenance.

The evaluated mask is the intersection of mapped, observed Cα positions in the
reference, AUTO primary and DEEPMSA primary. It is frozen before fitting and used
for both conditions, without score-driven pruning. `residue_mapping.csv` records
identifiers and exclusions; `frozen_mask.csv` records evaluated target positions.
At least three common positions are required. Coverage is reported against both
the full submitted sequence and the declared reference-normalization length.

TM-score and GDT-TS require an explicit `reference.reference_length`: the confirmed
full experimental-chain normalization length. It is never inferred from target
length or the observed coordinate count. Without it, RMSD and common-mask lDDT can
still run, while TM/GDT and reference-length coverage remain unavailable. The
choice is recorded, must be identical across conditions, and must be at least the
common-mask length. It does not imply every reference residue has coordinates.

## Structural accuracy

| Output | Definition and units | Missing-coordinate policy |
|---|---|---|
| `rmsd_angstrom` | Square root of the mean squared Cα distance after a NumPy SVD/Kabsch rigid fit. Rotation determinant is +1; reflections, outlier rejection and pruning are disallowed. Å; lower is better. | Same common mask for both conditions; divisor is mask size. Coverage is separate. |
| `lddt_ca_common_mask` | For unordered, distinct-residue reference Cα pairs closer than 15 Å, average the fractions whose absolute distance errors are strictly less than 0.5, 1, 2 and 4 Å. No superposition. Range 0–1; higher is better. | Reference contacts only within the frozen common mask. No eligible contacts gives missing score. The low-level function treats a missing mobile coordinate as an unconserved contact; the paired workflow excludes positions missing in either model from both masks. |
| `tm_score` | Optimized score from the official TMscore executable with fixed residue-index correspondence and `-l` declared reference length. The standard expression maximizes `sum_i 1/(1+(d_i/d0)^2) / L_ref`; official length-dependent d0 and search are used. Range 0–1; higher is better. | Only common-mask coordinates enter the numerator; the full declared normalization length remains in the denominator. No independent structural realignment of correspondences. |
| `gdt_ts` | Optimized GDT-TS returned by official TMscore: average counts within its 1, 2, 4 and 8 Å criteria, divided by declared reference length. Range 0–1, not percent; higher is better. | Official output divides by its input native-mask length even when `-l` is used. The wrapper multiplies that output by `mask_size/L_ref`, preserving full declared normalization. |
| `molprobity` | Separate all-atom stereochemical validation, currently unavailable. | No replacement score is generated. |

The reported lDDT is **Cα-only on a common mask**, not all-atom lDDT, whole-target
lDDT, or pLDDT. No stereochemical penalties or side-chain symmetry checks are
performed. Restricting the mask can raise apparent similarity; inspect coverage
and the excluded residue table. RMSD is likewise sensitive to evaluated coverage.
TM-score/GDT normalization keeps unscored positions in their denominator, but
coordinate absence and construct differences still limit interpretation.

## Optional standard-score tool

Set the confirmed full experimental-chain `reference.reference_length` and
provide an existing official **TMscore** executable from pylelab/USalign using
`tmscore_executable` in configuration, or place it on PATH as `TMscore`. This is
not the TMalign executable. The workflow does not download, compile or install it.
The wrapper runs a bounded CPU subprocess with argument-list invocation:

```text
TMscore AUTO_common_mask.pdb reference_common_mask.pdb -l L_ref
```

The derived PDB files have only Cα atoms, renumbered to identical target positions;
original structures are untouched. They store PDB coordinates to 0.001 Å, so
external scores can differ slightly from calculations on original CIF precision.
The wrapper requires the explicit user-normalized TM-score output, checks native
length and matched-residue count, captures stdout/stderr, and records executable
SHA-256. Missing binary, parser mismatch, timeout or tool failure gives unavailable
TM/GDT with a diagnostic. It never reports a single Kabsch-fit TM-like expression
or threshold fraction as optimized TM-score/GDT-TS.

The tests include synthetic documented-format parser tests. Live comparison with
the actual tool is separately skipped unless `PROJECT1_TMSCORE` points to an
existing executable; when enabled it compares wrapper output with direct official
tool output on temporary synthetic coordinates. A skip is a validation limitation,
not evidence that the external implementation was run. MolProbity/Phenix geometry
analysis has no automatic wrapper in this workflow and remains unavailable.

## Confidence and visualization

Public-server confidence is separate from experimental accuracy. pLDDT is 0–100,
pTM is 0–1, PAE is in Å, and ranking_score is only a ranking criterion. The import
code preserves CIF atom mapping and token chain/residue mapping. An atom mean is
reported distinctly from a residue mean: atoms are averaged within each residue,
then residue means are equally averaged. A mismatched/unavailable atom mapping
leaves the residue mean missing. No unavailable confidence field is imputed.

The Cα-error chart uses both models' independent rigid fits on the same fixed
mask. The ChimeraX script opens derived common-mask Cα structures and fits using
the saved target-position numbering without distance cutoffs. It is a mapped
Cα overlay, not an all-atom validation. Rendering requires a separately available
ChimeraX installation; missing rendering does not prevent tabular analysis.

Paired differences are `DEEPMSA minus AUTO`. A negative RMSD difference favors
DEEPMSA; positive TM/GDT/lDDT differences favor it. Server samples from one job are
not independent proteins. The three targets provide a descriptive comparison,
not general superiority, novelty or reconstruction of the blind CASP15 contest.
Disabling templates does not establish independence from training data. One
matched seed controls a setting but provides no replication.

## Sources and software provenance

Primary documentation was read on 2026-09-22:

- [Official TMscore command implementation](https://github.com/pylelab/USalign/blob/master/TMscore.cpp) and [score output/normalization](https://github.com/pylelab/USalign/blob/master/TMscore.h).
- [OpenStructure lDDT definitions](https://openstructure.org/docs/2.11/mol/alg/lddt/), including backbone-only Cα selection, radius, strict distance-error thresholds and contact denominator.
- [Public AlphaFold Server outputs](https://www.ebi.ac.uk/training/online/courses/alphafold/alphafold-3-and-alphafold-server/alphafold-server-your-gateway-to-alphafold-3/interpreting-results-from-alphafold-server/).
- [ChimeraX align command](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/align.html).
- Course materials audited in the parent workspace: sequence-alignment lecture pp.102–103 and pp.156–158; TBM lecture p.89; Nancy workflow Q6. No confirmed target/reference assignments were found in those sources.

NumPy supplies array operations/SVD; Biopython supplies BLOSUM62, sequence
alignment and coordinate readers; Matplotlib supplies standalone plots. Exact
tested package versions are captured in the review's dependency-version artifact
and test log, rather than inferred from the permitted ranges in requirements.txt.
