# Matthew's local analysis workflow

The current study uses **AlphaFold2 through ColabFold v1.6.3** on T1188, T1106s1 and T1122. The main comparison is automatic MMseqs2 versus custom JackHMMER for all three targets; T1188 also has a DeepMSA2 condition. Nancy's predictions, references, method notes and reported scores are the study inputs. This workflow adds reproducible local MSA analysis, coordinate comparisons, tables and figures.

Start with the [published results](results/matthew_analysis/REPORT.md) and [figure guide for Allen and Yu Wu](RESULTS_GUIDE.md). Nancy's [methods](docs/methods.md), [limitations](docs/limitations.md), [analysis summary](analysis/analysis_summary.md) and original files remain the source record.

## Run on Windows

Use Python 3.10 or newer. Open PowerShell in the repository root, then run:

```powershell
py -3 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m ipykernel install --user --name project1 --display-name "Project 1 (CPU)"
```

Open [Project1_Server_Workflow.ipynb](Project1_Server_Workflow.ipynb), select **Project 1 (CPU)**, then **Restart Kernel and Run All**. It is a code-only notebook. Despite its retained filename, its current configuration analyzes the uploaded AF2 study.

The equivalent command is:

```powershell
.\.venv\Scripts\python.exe -m project1 --mode uploaded
```

Run the software tests with:

```powershell
.\.venv\Scripts\python.exe -m pytest -q
```

The local workflow uses CPU analysis packages. It does not install or run AlphaFold or DeepMSA, submit website jobs, or require additional predictions. The tests use separate synthetic fixtures; study tables and figures use the real uploaded files.

## Files and outputs

The `uploaded_analysis` section of [config.json](config.json) maps the exact target sequences, seven alignments, seven prediction structures, three experimental references and Nancy's score table. It is the configuration for the current study.

- [results/matthew_analysis/REPORT.md](results/matthew_analysis/REPORT.md): stable entry point to the results.
- [combined_results.csv](results/matthew_analysis/combined_results.csv): condition-level table joining local analyses with team-reported scores.
- [paired_custom_minus_auto.csv](results/matthew_analysis/paired_custom_minus_auto.csv): paired differences for each target; negative RMSD differences favor custom, while positive lDDT/TM-score differences favor custom.
- [alignment_comparison.csv](results/matthew_analysis/alignment_comparison.csv): alignment depth, diversity, coverage and scoring summaries.
- `results/matthew_analysis/alignments/<target>/<condition>/`: alignment profiles, overview figures and supporting data.
- `results/matthew_analysis/structures/<target>/`: residue mapping and per-residue error/confidence figures.

Detailed run snapshots remain under `outputs/analyses/uploaded/`. The published `results/matthew_analysis/` folder provides stable paths for the team. Input files are read without overwriting Nancy's `analysis/`, `figures/`, `msa/`, `structures/`, `targets/`, notebooks or method notes.

## Reading the analysis

Local calculations describe the alignments and compare the supplied structures over the same mapped experimental residues within each target. Nancy's TM-score and pTM values are retained as team-reported results with their source. Independent reproduction of every reported metric is not a prerequisite for using the team's completed study.

Scored coverage is **573/630** residues for T1188, **71/122** for T1106s1 and **230/241** for T1122. Accuracy claims apply to those mapped regions. pLDDT and pTM describe confidence; RMSD, Cα-lDDT and TM-score describe agreement with experiment. These three targets support a descriptive case study.

The old AF3 submission preparation and its T1151s2 configuration are preserved as an earlier proposal. They are not instructions to start more jobs for the current study. See [analysis definitions](docs/UPLOADED_ANALYSIS.md) for calculation details.
