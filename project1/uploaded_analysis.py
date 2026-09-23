"""Reproduce analyses supported by Nancy's uploaded ColabFold study artifacts.

This is an explicitly separate evidence set from the earlier public-AF3-server
preparation policy. It never starts predictions or changes any uploaded file.
"""
from __future__ import annotations

from datetime import datetime, timezone
from importlib import metadata
import json
from pathlib import Path
import subprocess
import traceback
import zipfile

from .inputs import file_hash, read_fasta, sha256, write_json, write_text
from .workflow import load_config, resolve_path, write_csv


def uploads_available(config_path="config.json"):
    config = load_config(config_path)
    options = config.get("uploaded_analysis", {})
    return bool(options.get("enabled", True) and options.get("targets") and
                resolve_path(options.get("source_fasta", "targets/target_sequences.fasta"), config["_repo_root"]).is_file())


def ensure_uploaded_source(config_path="config.json"):
    """Do not silently turn an incomplete current checkout into the old study."""
    config = load_config(config_path)
    options = config.get("uploaded_analysis", {})
    if options.get("enabled", True) and options.get("targets"):
        source = resolve_path(options.get("source_fasta", "targets/target_sequences.fasta"), config["_repo_root"])
        if not source.is_file():
            raise FileNotFoundError(f"Configured AF2 study FASTA is missing: {source}. Restore the repository inputs before running; no legacy study was prepared.")


def _git(root, *args):
    completed = subprocess.run(["git", "-C", str(root), *args], capture_output=True, text=True, timeout=20)
    return completed.stdout.strip() if completed.returncode == 0 else "UNAVAILABLE"


def _inputs(config):
    root = config["_repo_root"]
    options = config["uploaded_analysis"]
    fasta = resolve_path(options["source_fasta"], root)
    records = read_fasta(fasta)
    targets, files = {}, [fasta]
    for ident, spec in options["targets"].items():
        matching = [row for row in records if row["id"] == spec["fasta_record_id"]]
        if len(matching) != 1 or matching[0]["errors"]:
            raise ValueError(f"{ident}: configured FASTA record absent, duplicated, or invalid")
        row = dict(spec, sequence=matching[0]["sequence"], source_header=matching[0]["header"])
        row["reference_path"] = str(resolve_path(spec["reference_path"], root))
        files.append(Path(row["reference_path"]))
        row["conditions"] = {}
        for condition, condition_spec in spec["conditions"].items():
            entry = dict(condition_spec)
            for key in ("alignment_path", "model_path"):
                entry[key] = str(resolve_path(entry[key], root))
                files.append(Path(entry[key]))
            row["conditions"][condition] = entry
        targets[ident] = row
    for item in ("targets/target_selection.md", "analysis/analysis_summary.md", "analysis/casp15_structure_comparison_metrics.csv", "docs/methods.md", "docs/limitations.md"):
        files.append(root / item)
    files.extend(sorted((root / "msa").glob("*/README.md")))
    inventory = [{"path": str(p), "relative_path": p.relative_to(root).as_posix(), "exists": p.is_file(),
                  "size_bytes": p.stat().st_size if p.is_file() else None,
                  "sha256": file_hash(p) if p.is_file() else None} for p in sorted(set(files))]
    return targets, inventory


def run_uploaded(config_path="config.json", package=True):
    config = load_config(config_path)
    root = config["_repo_root"]
    targets, inventory = _inputs(config)
    source_hashes = {p.name: file_hash(p) for p in sorted((root / "project1").glob("*.py"))}
    fingerprint = sha256(json.dumps({"inventory": inventory, "config": config["uploaded_analysis"], "code": source_hashes}, sort_keys=True).encode())[:16]
    output = root / "outputs" / "analyses" / "uploaded" / fingerprint
    output.mkdir(parents=True, exist_ok=True)
    result = {"status": "COMPLETED_FOR_AVAILABLE_UPLOADS", "mode": "uploaded_colabfold_evidence_audit",
              "source_git_revision": _git(root, "rev-parse", "HEAD"), "generated_utc": datetime.now(timezone.utc).isoformat(),
              "output_directory": str(output), "input_inventory": inventory, "code_sha256": source_hashes,
              "declared_method": "AlphaFold2 through ColabFold, automatic MMseqs2 versus custom HMMER; T1188 additionally DeepMSA2",
              "controls_status": ("TEAM_REPORTED_ACCEPTED: independent reconstruction is outside the agreed project scope" if config["uploaded_analysis"].get("team_results_accepted") else "UNVERIFIED: independent run reconstruction has not been performed"),
              "assignment_compatibility": ("AF2_ACCEPTED_BY_PROJECT_OWNER" if config["uploaded_analysis"].get("af2_accepted") else "NOT_RECORDED"),
              "casp_reference_assignments": {"status": "VERIFIED_OFFICIAL_TARGET_LIST", "source": "https://predictioncenter.org/casp15/targetlist.cgi?view=regular",
                                             "retrieved_date": "2026-09-22", "entries": {"T1188": "All groups/Ligand, A1, 630 residues, 8c6z", "T1106s1": "All groups, A1, 122 residues, 7qih", "T1122": "All groups, A1, 241 residues, 8bbt"}},
              "alignments": [], "issues": [], "targets": {}}
    write_json(output / "input_inventory.json", inventory)
    write_csv(output / "input_inventory.csv", inventory, ["relative_path", "path", "exists", "size_bytes", "sha256"])
    from .uploaded_msa import analyze_uploaded_msa
    for ident, target in targets.items():
        result["targets"][ident] = {"sequence_length": len(target["sequence"]), "source_header": target["source_header"],
                                     "reference_path": target["reference_path"], "reference_chain": target["reference_chain"]}
        for condition, spec in target["conditions"].items():
            print(f"Analyzing uploaded alignment: {ident} / {condition}")
            dest = output / "alignments" / ident / condition
            try:
                stats = analyze_uploaded_msa(Path(spec["alignment_path"]), target["sequence"], dest)
                stats.update({"target": ident, "condition": condition, "analysis_directory": str(dest)})
                result["alignments"].append(stats)
            except Exception as exc:
                issue = f"{ident}/{condition}: {type(exc).__name__}: {exc}"
                result["issues"].append(issue)
                write_text(output / "diagnostics" / f"{ident}_{condition}.txt", traceback.format_exc())
    summary_fields = ["target", "condition", "status", "query_length", "sequence_count_including_query", "unique_query_coordinate_rows",
                      "distinct_non_query_sequences", "duplicate_query_rows", "duplicate_fraction", "gap_fraction",
                      "mean_identity_to_query_excluding_query", "blosum62_normalized_sum_of_pairs", "blosum62_eligible_pair_denominator",
                      "neff", "neff_status", "source_sha256"]
    write_csv(output / "alignment_comparison.csv", result["alignments"], summary_fields)
    _plot_alignment_comparison(result["alignments"], output)
    from .uploaded_structures import audit_structures
    try:
        result["structures"] = audit_structures(root, output / "structures", targets)
        result["issues"].extend(result["structures"].get("issues", []))
    except Exception as exc:
        result["structures"] = {"status": "BLOCKED", "rows": [], "issues": [f"{type(exc).__name__}: {exc}"]}
        result["issues"].extend(result["structures"]["issues"])
        write_text(output / "diagnostics" / "structures.txt", traceback.format_exc())
    # Recheck original hashes; even a concurrent input change invalidates this snapshot.
    changed = [item["relative_path"] for item in inventory if item["exists"] and
               (not Path(item["path"]).is_file() or file_hash(item["path"]) != item["sha256"])]
    result["originals_integrity"] = "UNCHANGED" if not changed else "CHANGED_DURING_ANALYSIS"
    result["changed_inputs"] = changed
    if result["issues"] or changed:
        result["status"] = "PARTIAL_REQUIRES_REVIEW"
    versions = {"python": __import__("sys").version}
    for module in ("numpy", "biopython", "matplotlib", "nbformat", "nbclient", "pytest"):
        versions[module] = metadata.version(module)
    write_json(output / "provenance.json", {"git_revision": result["source_git_revision"], "git_status": _git(root, "status", "--short"),
                                            "versions": versions, "code_sha256": source_hashes, "config": config["uploaded_analysis"]})
    report = _report(result, output)
    write_text(output / "REPORT.md", report)
    write_json(output / "results.json", result)
    if config["uploaded_analysis"].get("publish_results", False):
        from .final_results import publish_results
        result["published_directory"] = str(publish_results(root, result, output))
        write_json(output / "results.json", result)
    write_json(root / "outputs" / "uploaded_analysis_current.json", {"output_directory": str(output), "report": str(output / "REPORT.md"),
                                                                    "source_git_revision": result["source_git_revision"], "fingerprint": fingerprint})
    write_text(root / "UPLOADED_RESULTS.md", "# Analysis of the current uploaded study\n\n" +
               f"Source revision: `{result['source_git_revision']}`. Original inputs: {result['originals_integrity']}.\n\n" +
               ("[Open the final project report](results/matthew_analysis/REPORT.md). [Slide and figure guide](RESULTS_GUIDE.md).\n\n" if result.get("published_directory") else f"[Open the current report]({(output / 'REPORT.md').relative_to(root).as_posix()}).\n\n") +
               "Run all cells in `Project1_Server_Workflow.ipynb` or `.\\.venv\\Scripts\\python.exe -m project1 --mode uploaded` to reproduce.\n\n" +
               "The final study uses Nancy's T1188, T1106s1 and T1122 sequences and AF2/ColabFold predictions. Team-reported results are accepted for this project. Earlier AF3 preparation files are historical.\n")
    if package:
        result["review_zip"] = str(_package(root, output, inventory))
    print(f"Uploaded-analysis report: {output / 'REPORT.md'}")
    print(f"Alignment conditions analyzed: {len(result['alignments'])}; original files: {result['originals_integrity']}")
    print("Study scope:", result["assignment_compatibility"], ";", result["controls_status"])
    if result.get("published_directory"):
        print("Final report and shareable figures:", result["published_directory"])
    return result


def _plot_alignment_comparison(rows, output):
    from matplotlib.figure import Figure
    if not rows:
        return
    labels = [f"{r['target']}\n{r['condition'].replace('_', ' ')}" for r in rows]
    fig = Figure(figsize=(12, 5.5), layout="constrained")
    ax = fig.subplots()
    x = list(range(len(rows)))
    ax.bar([i - .2 for i in x], [r["sequence_count_including_query"] for r in rows], width=.4, label="All records, including query")
    ax.bar([i + .2 for i in x], [r["unique_query_coordinate_rows"] for r in rows], width=.4, label="Distinct query-coordinate sequences")
    ax.set_yscale("log")
    ax.set_xticks(x, labels, rotation=25, ha="right")
    ax.set_ylabel("Sequence count (log scale)")
    ax.set_title("Uploaded alignments: raw depth and exact sequence diversity")
    ax.legend()
    for i, row in enumerate(rows):
        ax.text(i, max(row["sequence_count_including_query"], row["unique_query_coordinate_rows"]) * 1.12,
                f"{row['sequence_count_including_query']:,} / {row['unique_query_coordinate_rows']:,}", ha="center", fontsize=8)
    ax.set_ylim(.7, max(r["sequence_count_including_query"] for r in rows) * 3)
    fig.savefig(output / "alignment_depth_and_diversity.png", dpi=180)


def _report(result, output):
    rows = result["alignments"]
    lines = ["# Completed local analysis of the supplied AF2 study", "", f"Source Git revision: `{result['source_git_revision']}`.",
             f"Original-upload integrity after analysis: **{result['originals_integrity']}**.", "",
             "These are analyses of real uploaded MSAs and coordinates. Synthetic tests are separate. All colleague-uploaded files are preserved.", "",
             "## Newly generated alignment evidence", "", "| Target | Condition | Records incl. query | Distinct sequences | Distinct non-query sequences | State |",
             "|---|---|---:|---:|---:|---|"]
    for row in rows:
        lines.append(f"| {row['target']} | {row['condition']} | {row['sequence_count_including_query']:,} | {row['unique_query_coordinate_rows']:,} | {row.get('distinct_non_query_sequences', 'unavailable')} | {row['status']} |")
    structure_rows = result.get("structures", {}).get("rows", [])
    verified = sum(all(row.get(key) == "MATCH_WITHIN_ROUNDING" for key in
                       ("rmsd_report_comparison", "lddt_report_comparison", "plddt_report_comparison")) for row in structure_rows)
    coordinate_table = [f"**{verified}/{len(structure_rows)} coordinate comparisons reproduce all three reported RMSD, lDDT and mean pLDDT values within the uploaded table's rounding.**",
                        "", "| Target | Condition | Shared mask / target | RMSD (Å) | Cα lDDT (%) | Mean Cα pLDDT (whole model) |",
                        "|---|---|---:|---:|---:|---:|"]
    for row in structure_rows:
        coordinate_table.append(f"| {row['target']} | {row['condition']} | {row['scored_residue_count']}/{row['target_length']} | {row['rmsd_ca_angstrom']:.3f} | {100 * row['lddt_ca_common_mask']:.3f} | {row['ca_bfactor_residue_mean_all_mapped']:.3f} |")
    coordinate_table += ["", "The mask includes only positions mapped to the experimental structure and every compared model. In particular, T1106s1's comparison covers 71/122 residues; it does not establish whole-target accuracy. lDDT is a contact-weighted Cα calculation on this shared mask, not an all-atom score. Mean pLDDT uses all mapped model residues, not just the experimental mask.", ""]
    lines += ["", "See `alignment_comparison.csv`, `alignment_depth_and_diversity.png`, and `alignments/<target>/<condition>/` for exact coverage, gaps, identity, redundancy, normalized BLOSUM62 sum-of-pairs, conservation profiles, identity/coverage distributions, and labeled alignment overview figures.",
              "All-row statistics use complete alignments. The overview image displays a bounded, deterministic subset of rows only; it is not an analysis subsample. Exact Neff is unavailable above the declared 500-row bound. Reported DeepMSA Nf is a distinct service-supplied statistic and is not substituted for local Neff.",
              "T1122's three automatic records are identical copies of its query; its custom HMMER is query-only. A sequence-only outcome can be analyzed as a documented low-information condition. It is not treated as a homolog-rich alignment or used to generate an artificial DeepMSA request.", "",
              "## Coordinate and confidence audit", "", *coordinate_table,
              "`structures/` contains residue mapping/mask evidence, numerical checks against the uploaded table, and per-residue figures where coordinate identity could be established. All uploaded conditions for a target use one shared coordinate mask. Computed similarity remains distinct from verification of experimental controls.",
              "RMSD and common-mask C-alpha lDDT are independently recalculated from uploaded coordinates. Differences from colleague-reported values are retained. Uploaded TM-score and pTM are not independently verified by this calculation. ColabFold PDB B-values can support a declared pLDDT interpretation, with the aggregation convention stated; PAE matrices, pTM values and complete ranking cannot be recovered from coordinates alone.", "",
              "## Accepted project scope", "",
              "The final project uses T1188, T1106s1 and T1122, with automatic MMseqs2 versus custom HMMER as the main paired comparison and DeepMSA2 additionally for T1188.",
              f"Scope decision: {result['assignment_compatibility']}. Team-result policy: {result['controls_status']}.",
              "Nancy's supplied TM-score, pTM and other result-table values are retained as team-reported results. Independently reconstructing every prediction setting or obtaining raw confidence/search exports is not required to complete this agreed scope. The documented TM-scores use mapped-residue normalization.",
              "The [official CASP15 regular-target list](https://predictioncenter.org/casp15/targetlist.cgi?view=regular), checked 2026-09-22, confirms All-groups/A1 entries and PDB assignments for T1188/8C6Z, T1106s1/7QIH and T1122/8BBT. Coordinate mapping and experimental coverage are checked separately. T1106s1's membership in a biological complex does not invalidate its individual A1 CASP entry.", "",
              "## Interpretation limits", "",
              "This is a descriptive three-target case study. Search method, filtering, model choice and target biology may all matter; these results do not establish a universal or isolated causal effect of MSA depth. Optional reconstruction checks are not project-completion blockers.", "",
              "## Reproduction", "", "Run `Project1_Server_Workflow.ipynb` (all code cells), or `python -m project1 --mode uploaded`. Edit only `config.json` for explicit input/condition mappings. No server submissions, local predictions, downloads of scoring executables, or model training are performed.", "",
              "Software tests: `python -m pytest -q`. Test sequences and coordinates are isolated synthetic fixtures. `input_inventory.csv` records all original file sizes and hashes; `provenance.json` records software and source revision."]
    issues = result.get("issues", [])
    if issues:
        lines += ["", "## Issues requiring review", ""] + [f"- {issue}" for issue in issues]
    return "\n".join(lines) + "\n"


def _package(root, output, inventory):
    """Allowlisted snapshot; no Git internals, environments, older analyses or secrets."""
    from .review import _artifact_copy, redact_text
    dest = root / "outputs" / "reviews"
    dest.mkdir(parents=True, exist_ok=True)
    path = dest / ("uploaded_analysis_review_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ") + ".zip")
    selected = set(output.rglob("*"))
    selected.update((root / "project1").glob("*.py"))
    selected.update((root / "tests").glob("*.py"))
    selected.update((root / "results" / "matthew_analysis").rglob("*"))
    selected.update(root / item for item in ("Project1_Server_Workflow.ipynb", "config.json", "requirements.txt", "UPLOADED_RESULTS.md", "README.md", "MATTHEW_WORKFLOW.md", "RESULTS_GUIDE.md", "START_HERE.md", "docs/METRICS.md", "docs/UPLOADED_ANALYSIS.md", "docs/legacy_af3_config.json", "docs/VALIDATION.md", "outputs/tests.log"))
    selected.update(Path(item["path"]) for item in inventory if item["exists"])
    originals = {Path(item["path"]) for item in inventory if item["exists"]}
    scientific_suffixes = {".a3m", ".fasta", ".fa", ".faa", ".pdb", ".cif", ".mmcif"}
    excluded, entries, total = [], [], 0
    with zipfile.ZipFile(path, "x", compression=zipfile.ZIP_DEFLATED) as archive:
        for source in sorted(selected):
            if not source.is_file() or source.is_symlink():
                continue
            size = source.stat().st_size
            if size > 20 * 1024 * 1024 or total + size > 180 * 1024 * 1024:
                excluded.append({"path": str(source), "size_bytes": size, "sha256": file_hash(source), "reason": "Size bound; original remains intact, not truncated"})
                continue
            original_bytes = source.read_bytes()
            if source in originals and source.suffix.lower() in scientific_suffixes:
                # Preserve scientific bytes (including line endings) unless an
                # access-bearing URL genuinely needs redaction in a package copy.
                decoded = original_bytes.decode("utf-8")
                sanitized = redact_text(decoded)
                content = original_bytes if sanitized == decoded else sanitized.encode("utf-8")
            else:
                content = _artifact_copy(source)
            relative = source.relative_to(root).as_posix()
            archive.writestr(relative, content)
            entries.append({"path": relative, "original_sha256": sha256(original_bytes),
                            "packaged_sha256": sha256(content), "original_input": source in originals,
                            "copy_mode": "BYTE_IDENTICAL" if content == original_bytes else "SANITIZED_PACKAGE_COPY"})
            total += len(content)
        archive.writestr("EXCLUDED_ARTIFACTS.json", json.dumps(excluded, indent=2))
        archive.writestr("PACKAGE_MANIFEST.json", json.dumps(entries, indent=2))
    write_json(root / "outputs" / "uploaded_review_current.json", {"path": str(path), "sha256": file_hash(path)})
    return path
