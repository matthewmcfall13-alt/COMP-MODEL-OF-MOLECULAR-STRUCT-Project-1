# Start here — completed AF2 study

The current study uses **T1188, T1106s1 and T1122**, predicted with AlphaFold2 through ColabFold v1.6.3. The main comparison is automatic MMseqs2 versus custom JackHMMER; T1188 also has a DeepMSA2 condition. The team's alignments, structures and reported results are already in the repository. No additional server jobs are needed for this study.

1. Read the [project overview](README.md) and [current results report](results/matthew_analysis/REPORT.md).
2. Use the [results and figure guide](RESULTS_GUIDE.md) to locate tables, alignment figures and structural overlays.
3. To reproduce the local analysis, follow [MATTHEW_WORKFLOW.md](MATTHEW_WORKFLOW.md). Run **Restart Kernel and Run All** in the code-only [Project1_Server_Workflow.ipynb](Project1_Server_Workflow.ipynb), or run `python -m project1 --mode uploaded` in the configured Python environment.

Nancy's [methods](docs/methods.md), [limitations](docs/limitations.md) and [analysis summary](analysis/analysis_summary.md) document the completed experiment. [docs/UPLOADED_ANALYSIS.md](docs/UPLOADED_ANALYSIS.md) explains Matthew's added local calculations and the handling of team-reported metrics.

The older AF3 submission files and the T1151s2 pilot configuration belong to an earlier proposal. They are preserved for history and are not the next steps for this completed AF2 study.
