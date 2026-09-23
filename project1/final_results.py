"""Publish Matthew's derived results beside, without editing, team originals."""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
from pathlib import Path


OWNER = "project1.final_results.v1"
MARKER = ".generated_manifest.json"
CONDITIONS = {
    "T1188": ("Automatic MMseqs2", "Custom HMMER", "Custom DeepMSA2"),
    "T1106s1": ("Automatic MMseqs2", "Custom HMMER"),
    "T1122": ("Automatic MMseqs2", "Custom HMMER"),
}
TEAM_METRICS = {
    "MSA_depth": "team_msa_depth", "pLDDT": "team_plddt", "pTM": "team_ptm",
    "TM_score": "team_tm_score", "lDDT_CA": "team_lddt_ca_percent",
    "RMSD_CA_A": "team_rmsd_ca_angstrom", "RMSD100_A": "team_rmsd100_angstrom",
    "Mapped_residues": "team_mapped_residues", "Resolution_A": "resolution_angstrom",
}
MSA_METRICS = (
    "sequence_count_including_query", "unique_query_coordinate_rows", "distinct_non_query_sequences",
    "duplicate_query_rows", "duplicate_fraction", "gap_fraction", "mean_identity_to_query_excluding_query",
    "blosum62_normalized_sum_of_pairs", "neff", "neff_status",
)
COORDINATE_METRICS = (
    "target_length", "scored_residue_count", "scored_fraction_of_target", "rmsd_ca_angstrom",
    "lddt_ca_common_mask", "ca_bfactor_residue_mean_all_mapped",
)


def _digest(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _number(value):
    if value is None or value == "":
        return None
    try:
        number = float(value)
        return number if math.isfinite(number) else None
    except (TypeError, ValueError):
        return None


def _csv_bytes(rows: list[dict], fields=None) -> bytes:
    fields = fields or list(dict.fromkeys(key for row in rows for key in row))
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def _condition(name: str) -> str:
    return str(name).replace("_", " ")


def _key(row: dict, team=False):
    return (row.get("Target" if team else "target"),
            _condition(row.get("Condition" if team else "condition", "")))


def _index(rows, team=False):
    indexed = {}
    for row in rows:
        key = _key(row, team)
        if key in indexed:
            raise ValueError(f"Duplicate target/condition cannot be merged: {key}")
        indexed[key] = row
    return indexed


def _merge(result: dict, team_rows: list[dict]) -> list[dict]:
    msas = _index(result.get("alignments", []))
    structures = _index(result.get("structures", {}).get("rows", []))
    team = _index(team_rows, team=True)
    expected = [(target, condition) for target, conditions in CONDITIONS.items() for condition in conditions]
    extras = sorted((set(msas) | set(structures) | set(team)) - set(expected))
    merged = []
    for target, condition in expected + extras:
        key = (target, condition)
        alignment, coordinates, reported = msas.get(key, {}), structures.get(key, {}), team.get(key, {})
        missing = [name for name, data in (("ALIGNMENT", alignment), ("COORDINATES", coordinates), ("TEAM_RESULTS", reported)) if not data]
        row = {"target": target, "condition": condition,
               "availability": "COMPLETE" if not missing else "MISSING_" + "+".join(missing),
               "difficulty": reported.get("Difficulty"), "experimental_pdb": reported.get("PDB"),
               "team_selected_model": reported.get("Rank1_model"),
               "alignment_status": alignment.get("status")}
        row.update({name: alignment.get(name) for name in MSA_METRICS})
        row.update({destination: _number(reported.get(source)) for source, destination in TEAM_METRICS.items()})
        row.update({"calculated_" + name: coordinates.get(name) for name in COORDINATE_METRICS})
        row["team_results_source"] = "analysis/casp15_structure_comparison_metrics.csv" if reported else None
        row["alignment_source_sha256"] = alignment.get("source_sha256")
        merged.append(row)
    return merged


def _pairs(rows):
    indexed = {(r["target"], r["condition"]): r for r in rows}
    deltas = []
    fields = list(TEAM_METRICS.values()) + [name for name in MSA_METRICS if name != "neff_status"]
    fields += ["calculated_" + name for name in COORDINATE_METRICS]
    for row in rows:
        if row["condition"] == "Automatic MMseqs2":
            continue
        auto = indexed.get((row["target"], "Automatic MMseqs2"), {})
        paired = {"target": row["target"], "custom_condition": row["condition"],
                  "comparison": "Custom minus Automatic MMseqs2",
                  "availability": "COMPLETE" if row["availability"] == auto.get("availability") == "COMPLETE" else "PARTIAL",
                  "scope": "Main comparison" if row["condition"] == "Custom HMMER" else "Additional comparison"}
        for name in fields:
            custom_value, auto_value = _number(row.get(name)), _number(auto.get(name))
            paired[name + "_delta"] = custom_value - auto_value if custom_value is not None and auto_value is not None else None
        deltas.append(paired)
    return deltas


def _portable_value(value: str, root: Path, output: Path) -> str:
    # Derived CSVs can refer to the original hash directory. Replace that prefix
    # before the repository prefix so paths stay relative to this publication.
    original = value
    for source, replacement in ((output, "."), (root, "../..")):
        for spelling in (str(source), source.as_posix()):
            value = value.replace(spelling + "\\", replacement + "/").replace(spelling + "/", replacement + "/")
            if value == spelling:
                value = replacement
    return value.replace("\\", "/") if value != original else value


def _inventory(root: Path, result: dict):
    rows = []
    seen = set()
    for item in result.get("input_inventory", []):
        path = Path(item["path"]).resolve()
        if not path.is_relative_to(root):
            raise ValueError("Published source inventory requires repository-contained inputs")
        relative = path.relative_to(root).as_posix()
        seen.add(relative)
        rows.append({"path": relative, "exists": bool(item.get("exists")),
                     "size_bytes": item.get("size_bytes"), "sha256": item.get("sha256")})
    # These newly uploaded notes establish the accepted team's actual method.
    for relative in ("analysis/casp15_structure_comparison_metrics.csv", "docs/methods.md", "docs/limitations.md"):
        path = root / relative
        if relative not in seen and path.is_file():
            content = path.read_bytes()
            rows.append({"path": relative, "exists": True, "size_bytes": len(content), "sha256": _digest(content)})
    return sorted(rows, key=lambda item: item["path"])


def _format(value, precision=2):
    number = _number(value)
    return "—" if number is None else f"{number:.{precision}f}"


def _findings(rows, pairs):
    lines = []
    for pair in pairs:
        if pair["availability"] != "COMPLETE":
            continue
        rmsd, tm, lddt = (pair.get(name) for name in ("team_rmsd_ca_angstrom_delta", "team_tm_score_delta", "team_lddt_ca_percent_delta"))
        lines.append(f"- **{pair['target']}, {pair['custom_condition']} versus automatic:** reported changes are {_format(tm, 3)} TM-score, {_format(lddt)} Cα-lDDT percentage points, and {_format(rmsd)} Å RMSD (custom minus automatic). Higher TM/lDDT and lower RMSD indicate better agreement on the stated mask.")
    hard = [row for row in rows if row["target"] == "T1122" and row.get("distinct_non_query_sequences") is not None]
    if hard and all(row["distinct_non_query_sequences"] == 0 for row in hard):
        lines.append("- **T1122 has no distinct sequence beyond its query in the supplied alignments.** Three automatic records are duplicate query sequences; the HMMER alignment is query-only. Their counts do not provide three independent homologs.")
    moderate = next((row for row in rows if row["target"] == "T1106s1" and row.get("calculated_scored_residue_count") is not None), None)
    if moderate:
        lines.append(f"- **T1106s1 coverage:** {moderate['calculated_scored_residue_count']}/{moderate['calculated_target_length']} target residues are scored. Its structural results describe the experimentally mapped segment, not the entire target.")
    return lines


def _report(result, rows, pairs, complete):
    completed = sum(row["availability"] == "COMPLETE" for row in rows)
    lines = ["# Matthew's alignment and structure results", "",
             "**Status: " + ("complete for the seven supplied study conditions." if complete else "partial; available results are retained and missing inputs are listed below.") + "**", "",
             f"Source Git revision: `{result.get('source_git_revision', 'not recorded')}`. Generated: {result.get('generated_utc', 'not recorded')}. {completed}/{len(rows)} condition rows combine an analyzed alignment, mapped coordinates and the team's result table.", "",
             "The accepted study uses Nancy's T1188, T1106s1 and T1122 targets with AlphaFold2/ColabFold. The main comparison is automatic MMseqs2 versus custom HMMER; T1188 additionally compares DeepMSA2. Nancy's structural results are accepted team results and retained in the `team_` columns. Matthew's alignment statistics and coordinate recalculations appear separately; neither replaces the original measurements.", "",
             "[Team methods](../../docs/methods.md) · [Team limitations](../../docs/limitations.md) · [Original result table](../../analysis/casp15_structure_comparison_metrics.csv) · [Original structure figures](../../figures/)", "",
             "## Combined results", "",
             "| Target | Condition | Rows / distinct sequences | TM-score | Cα-lDDT (%) | RMSD (Å) | Available inputs |",
             "|---|---|---:|---:|---:|---:|---|"]
    for row in rows:
        count = " / ".join(_format(row.get(field), 0) for field in ("sequence_count_including_query", "unique_query_coordinate_rows"))
        lines.append(f"| {row['target']} | {row['condition']} | {count} | {_format(row['team_tm_score'], 3)} | {_format(row['team_lddt_ca_percent'])} | {_format(row['team_rmsd_ca_angstrom'])} | {row['availability']} |")
    lines += ["", "Structural values in this table are Nancy's accepted reported values; complete precision for recalculated values is in `combined_results.csv`. Blank values remain missing, never zero.", "",
              "## What the comparison shows", ""] + _findings(rows, pairs)
    lines += ["", "These are descriptive within-target comparisons. Search method, filtering, redundancy, coverage and model selection vary alongside depth; this three-target case study does not isolate a causal effect of sequence count. Small numerical differences are not evidence of statistical significance.", "",
              "## Files and metric conventions", "",
              "- `combined_results.csv`: joined alignment statistics, accepted team structural/confidence metrics and separate coordinate recalculations for every planned condition.",
              "- `paired_custom_minus_auto.csv`: signed differences. Negative RMSD/RMSD100 and positive TM-score/lDDT favor custom; pLDDT/pTM measure confidence, not experimental accuracy.",
              "- `source_inventory.csv`: repository-relative source paths and original SHA-256 hashes. Team inputs are read only.",
              "- `alignment_depth_and_diversity.png` and `alignments/`: full-row coverage, gaps, query identity, redundancy, BLOSUM62 sum-of-pairs, entropy, and labeled alignment overviews. Overview rows are a deterministic display subset; numerical statistics use all rows.",
              "- `structures/`: shared residue masks, residue mappings, per-residue error/confidence/coverage figures and recalculated comparisons.",
              "- `FIGURE_CAPTIONS.md`: text Allen can use with these figures independently of a slide deck.", "",
              "Counts include the query. Distinct sequences are exact insertion-stripped query-coordinate strings, including gap positions; they are not a count of independent evolutionary observations. Entropy excludes gaps/noncanonical symbols. BLOSUM62 sum-of-pairs excludes pairs containing gaps/noncanonical residues and is normalized by eligible pair count. Exact local Neff is only calculated within the documented 500-row limit; a blank value for a large alignment is unavailable, not zero, and is distinct from DeepMSA's reported Nf.", "",
              "Coordinate recalculations use the same mapped Cα mask for all conditions of each target. Kabsch RMSD uses no outlier pruning. Cα-lDDT uses reference contacts under 15 Å and strict distance-error thresholds 0.5/1/2/4 Å. Whole-model mean pLDDT uses one value per mapped model residue. Nancy's TM-score normalization is the mapped experimental count; confidence and accuracy remain separate.", "",
              "## Reproduction", "", "From the repository, run `python -m project1 --mode uploaded` or all code cells in `Project1_Server_Workflow.ipynb`. Publication is enabled by `uploaded_analysis.publish_results` in `config.json`. Only files owned by this generated folder are refreshed; original team folders are preserved."]
    missing = [f"- {row['target']} / {row['condition']}: {row['availability']}" for row in rows if row["availability"] != "COMPLETE"]
    if missing or result.get("issues") or result.get("originals_integrity") != "UNCHANGED":
        lines += ["", "## Availability and run notes", ""] + missing
        lines += ["- " + str(issue) for issue in result.get("issues", [])]
        lines.append("- Original-input integrity state: " + str(result.get("originals_integrity", "not recorded")))
    return "\n".join(lines) + "\n"


def _captions(rows, pairs):
    lines = ["# Figure captions and result notes for Allen", "",
             "**Alignment depth and diversity (`alignment_depth_and_diversity.png`).** All supplied records versus distinct query-coordinate sequence strings for each MSA condition. Counts include the query; exact duplicates retain their multiplicity in the full-row statistics. The logarithmic axis displays the large depth range without treating sequence count as a direct accuracy measure.", "",
             "**Alignment profiles (`alignments/<target>/<condition>/alignment_profiles.png`).** Per-position non-gap coverage and gap fraction use all records. Canonical-residue entropy and consensus exclude gaps/ambiguous residues; non-query profiles exclude row zero. Identity/coverage distributions describe all eligible rows after the query. Query-only or duplicate-query conditions contain no distinct homolog-sequence information.", "",
             "**Alignment overview (`alignments/<target>/<condition>/alignment_overview.png`).** Query-coordinate view of up to 80 deterministically spaced file-order rows. Colors distinguish gaps, noncanonical residues, mismatches and query identity. This is a display subset; every input row contributes to numerical summaries.", "",
             "**Structure profiles (`structures/<target>/per_residue_error_confidence_coverage.png`).** Cα deviations after independent rigid fits on one shared experimental mask, predicted per-residue pLDDT, and coordinate coverage. Missing experimental positions are unscored. pLDDT describes model confidence rather than measured structural accuracy.", "",
             "## Target-specific notes", ""] + _findings(rows, pairs)
    lines += ["", "Nancy's existing overview/accuracy figures remain in `../../figures/`; cite `../../analysis/casp15_structure_comparison_metrics.csv` for their reported values. The new figures supplement those results with alignment composition and residue-level evidence. Main comparisons are HMMER versus automatic; the T1188 DeepMSA2 result is an additional comparison."]
    return "\n".join(lines) + "\n"


def publish_results(root, result: dict, output) -> Path:
    """Write portable derived artifacts to ``results/matthew_analysis``.

    The returned Path is the stable directory. Original files are read only.
    An ownership manifest prevents replacement of unrelated or edited files;
    stale files are removed only when their bytes still match our last write.
    """
    root, output = Path(root).resolve(), Path(output).resolve()
    destination = root / "results" / "matthew_analysis"
    if destination.is_symlink() or (root / "results").is_symlink():
        raise ValueError("Publication destination must not be a symlink")
    if output == destination or not output.is_relative_to(root) or not output.is_dir():
        raise ValueError("Require an existing separate derived-output folder inside the repository")
    previous = {}
    marker = destination / MARKER
    if destination.exists():
        if marker.is_file():
            manifest = json.loads(marker.read_text(encoding="utf-8"))
            if manifest.get("owner") != OWNER:
                raise ValueError("Publication folder has a different owner")
            previous = manifest.get("files", {})
        elif any(destination.iterdir()):
            raise ValueError("Refusing to overwrite a nonempty publication folder without its ownership manifest")
    for relative, digest in previous.items():
        owned = (destination / relative).resolve()
        if not owned.is_relative_to(destination.resolve()) or owned.is_symlink():
            raise ValueError("Unsafe path in publication manifest")
        if owned.is_file() and _digest(owned.read_bytes()) != digest:
            raise ValueError(f"Generated file was edited; preserve or move it before refreshing: {relative}")

    table_path = root / "analysis" / "casp15_structure_comparison_metrics.csv"
    with table_path.open(encoding="utf-8-sig", newline="") as handle:
        team_rows = list(csv.DictReader(handle))
    rows = _merge(result, team_rows)
    pairs = _pairs(rows)
    complete = (len(rows) == 7 and all(row["availability"] == "COMPLETE" for row in rows)
                and result.get("originals_integrity") == "UNCHANGED"
                and not result.get("issues") and not str(result.get("status", "")).startswith(("PARTIAL", "BLOCKED")))
    artifacts = {}
    for source in sorted(output.rglob("*")):
        if not source.is_file() or source.is_symlink() or source.suffix.lower() not in {".csv", ".png"}:
            continue
        relative = source.relative_to(output).as_posix()
        if relative == "input_inventory.csv":
            continue
        if source.suffix.lower() == ".png":
            artifacts[relative] = source.read_bytes()
        else:
            with source.open(encoding="utf-8-sig", newline="") as handle:
                reader = csv.DictReader(handle)
                fields = reader.fieldnames or []
                portable = [{key: _portable_value(value or "", root, output) for key, value in row.items()} for row in reader]
            artifacts[relative] = _csv_bytes(portable, fields)
    artifacts["combined_results.csv"] = _csv_bytes(rows)
    artifacts["paired_custom_minus_auto.csv"] = _csv_bytes(pairs)
    artifacts["source_inventory.csv"] = _csv_bytes(_inventory(root, result), ["path", "exists", "size_bytes", "sha256"])
    artifacts["REPORT.md"] = _report(result, rows, pairs, complete).encode("utf-8")
    artifacts["FIGURE_CAPTIONS.md"] = _captions(rows, pairs).encode("utf-8")
    for relative in artifacts:
        path = destination / relative
        if path.exists() and relative not in previous:
            raise ValueError(f"Refusing to overwrite an unowned publication file: {relative}")
        if not path.resolve().is_relative_to(destination.resolve()) or path.is_symlink():
            raise ValueError("Unsafe publication output path")
    destination.mkdir(parents=True, exist_ok=True)
    for relative, content in artifacts.items():
        path = destination / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
    for relative in previous.keys() - artifacts.keys():
        path = destination / relative
        if path.is_file():
            path.unlink()
    marker.write_text(json.dumps({"owner": OWNER, "source_git_revision": result.get("source_git_revision"),
                                 "complete": complete, "files": {name: _digest(content) for name, content in sorted(artifacts.items())}}, indent=2) + "\n", encoding="utf-8")
    return destination
