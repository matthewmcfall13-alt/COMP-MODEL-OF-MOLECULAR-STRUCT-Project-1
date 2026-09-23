# Results guide for Allen and Yu Wu

Use the [current results report](results/matthew_analysis/REPORT.md), [combined table](results/matthew_analysis/combined_results.csv) and [paired differences](results/matthew_analysis/paired_custom_minus_auto.csv). The completed study uses **AF2/ColabFold v1.6.3**: automatic MMseqs2 versus custom JackHMMER for T1188, T1106s1 and T1122, plus a DeepMSA2 comparison for T1188.

The existing [team analysis summary](analysis/analysis_summary.md), [methods](docs/methods.md) and [limitations](docs/limitations.md) provide the scientific narrative. These results do not require new prediction jobs or independent reruns of the team's reported scores.

## Figures to use

All paths below are relative to this repository. New figures are in `results/matthew_analysis/`; Nancy's original figures remain in `figures/`.

| Purpose | File | Caption or talking point |
|---|---|---|
| Compare raw MSA depth | [Original depth chart](figures/figure_msa_depth.png) | Sequence count varies widely; count alone does not measure independent evolutionary information. |
| Show alignment quality and coverage | `results/matthew_analysis/alignments/<target>/<condition>/alignment_profiles.png` | Coverage, gaps, conservation and identity describe different aspects of the alignment. |
| Show the actual alignment pattern | `results/matthew_analysis/alignments/<target>/<condition>/alignment_overview.png` | The labeled row subset is for display; reported statistics use the full alignment. |
| Compare global structural accuracy | [Original TM-score chart](figures/figure_tm_score.png) | Use the team's reported TM-scores, normalized to mapped experimental residues as described in its analysis summary. Higher is better. |
| Compare local structural accuracy | [Original Cα-lDDT chart](figures/figure_lddt_ca.png) | Local agreement can improve even when global RMSD worsens. Higher is better. |
| Compare size-normalized deviation | [Original RMSD100 chart](figures/figure_rmsd100.png) | Label this RMSD100, not raw RMSD. Lower is better. |
| Explain regional differences | `results/matthew_analysis/structures/<target>/per_residue_error_confidence_coverage.png` | Relate residue-level deviation to confidence and experimental coverage. |
| Easy target overlay | [T1188 overlay](figures/overlay_T1188_easy.png) | All three strategies closely match the resolved experimental structure. |
| Moderate target overlay | [T1106s1 overlay](figures/overlay_T1106s1_moderate.png) | The metric-dependent comparison covers only 71 of 122 target residues. |
| Hard target overlay | [T1122 overlay](figures/overlay_T1122_hard.png) | Both predictions remain substantially less accurate than the easy-target predictions. |
| Easy-target bonus comparison | [T1188 three-way TM-score](figures/figure_easy_three_way_tm_score.png) | The much deeper DeepMSA2 alignment produces only a small change for this target. |

Substitute `T1188`, `T1106s1` or `T1122` for `<target>`. The condition folders are `Automatic_MMseqs2` and `Custom_HMMER`, with `Custom_DeepMSA2` additionally available for T1188.

## Three findings to explain

1. **T1188: similar accuracy despite very different depths.** Automatic, HMMER and DeepMSA2 contain 11,337, 213 and 25,136 records respectively. The filtered HMMER alignment gives a similarly accurate structure; this is consistent with diminishing returns for this target.
2. **T1106s1: the answer depends on the metric.** Custom HMMER improves local Cα-lDDT while automatic MMseqs2 has lower global RMSD and higher reported TM-score. Explain both results and the limited experimental coverage instead of naming one universal winner.
3. **T1122: raw depth overstates the available information.** The three automatic records are identical copies of the query; the custom HMMER file contains only the query. Neither alignment has non-query sequence diversity. The custom model's modest improvement does not demonstrate that adding homologs caused an improvement.

## Labels to retain

- The main custom condition is **JackHMMER**, not DeepMSA2. DeepMSA2 is the T1188 bonus comparison.
- Custom JackHMMER used one iteration against UniProt 2025_01. T1188 hits were filtered to at least 70% target coverage; both alignment size and composition change with that filter.
- State mapped experimental coverage: T1188 **573/630**, T1106s1 **71/122**, T1122 **230/241**. These are not whole-target accuracy estimates.
- Distinguish confidence (**pLDDT, pTM**) from experimental accuracy (**TM-score, Cα-lDDT, RMSD/RMSD100**).
- Paired differences are **custom minus automatic**. Their favorable direction depends on the metric.
- Call the project a **three-target case study**. The results do not establish a universal causal effect of MSA depth.

For rerunning the analysis or using the code-only notebook, see [Matthew's workflow](MATTHEW_WORKFLOW.md). No additional server submission is needed.
