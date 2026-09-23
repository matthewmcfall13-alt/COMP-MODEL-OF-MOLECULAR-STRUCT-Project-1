"""Reproducible coordinate audit of uploaded models, without prediction services.

This audits the files actually supplied. Matching experimental controls, rank
selection, pTM and optimized TM-score cannot be established from a PDB alone.
All artifacts are derived copies; source structures and summaries are read only.
"""
from __future__ import annotations

import csv
from decimal import Decimal, InvalidOperation
import hashlib
import json
from pathlib import Path
import re

import numpy as np

from .metrics import (
    AA, MappingError, kabsch_rmsd, lddt_ca, map_to_target,
    metric_definitions, read_ca_structure,
)


CASP15_TARGET_LIST = "https://predictioncenter.org/casp15/targetlist.cgi?view=regular"
CASP15_ASSIGNMENTS = {"T1188": "8C6Z", "T1106s1": "7QIH", "T1122": "8BBT"}


def _hash(path: Path) -> str:
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest() if hasattr(hashlib, "file_digest") else hashlib.sha256(handle.read()).hexdigest()


def _csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    fields = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _path(root: Path, value) -> Path:
    path = Path(value)
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def rounding_comparison(value: float | None, reported: str | None) -> dict:
    """Compare to the precision actually printed, inclusive of rounding bounds."""
    if value is None or reported in (None, ""):
        return {"status": "UNAVAILABLE", "reported": reported, "difference": None, "tolerance": None}
    try:
        number = Decimal(str(reported))
        if not number.is_finite() or not np.isfinite(value):
            raise ValueError("Nonfinite score")
        tolerance = float(Decimal("0.5") * Decimal(10) ** number.as_tuple().exponent)
        difference = float(value) - float(number)
        return {"status": "MATCH_WITHIN_ROUNDING" if abs(difference) <= tolerance + 1e-10 else "DIFFERS",
                "reported": float(number), "difference": difference, "tolerance": tolerance}
    except (InvalidOperation, ValueError, OverflowError):
        return {"status": "INVALID_REPORTED_VALUE", "reported": reported, "difference": None, "tolerance": None}


def _header(path: Path) -> dict:
    """Local PDB header evidence; not a fresh external database verification."""
    text = path.read_text(encoding="utf-8", errors="replace")
    methods = [line[10:].strip() for line in text.splitlines() if line.startswith("EXPDTA")]
    resolution = re.search(r"REMARK\s+2 RESOLUTION\.\s+([0-9.]+) ANGSTROMS", text)
    header = next((line for line in text.splitlines() if line.startswith("HEADER")), "")
    return {"evidence_source": "Uploaded PDB header", "pdb_id": header[62:66].strip() or None,
            "experimental_method": " ".join(methods) or None,
            "resolution_angstrom": float(resolution.group(1)) if resolution else None}


def _condition_label(name: str, condition: dict) -> str:
    if condition.get("reported_condition"):
        return condition["reported_condition"]
    simple = name.lower().replace("_", "").replace(" ", "")
    return {"auto": "Automatic MMseqs2", "automaticmmseqs2": "Automatic MMseqs2",
            "hmmer": "Custom HMMER", "customhmmer": "Custom HMMER",
            "deepmsa": "Custom DeepMSA2", "deepmsa2": "Custom DeepMSA2",
            "customdeepmsa2": "Custom DeepMSA2"}.get(simple, name)


def _plots(dest: Path, target: str, length: int, conditions: list[str], mappings: dict,
           positions: list[int], profiles: dict[str, list[dict]]) -> str | None:
    try:
        from matplotlib.figure import Figure
    except ImportError:
        return None
    figure = Figure(figsize=(11, 8), layout="constrained")
    error_ax, confidence_ax, coverage_ax = figure.subplots(3, 1, gridspec_kw={"height_ratios": [2, 2, 1]})
    for condition in conditions:
        rows = profiles[condition]
        # NaN breaks lines at absent positions; no visual bridging of missing residues.
        error_ax.plot([r["target_position"] for r in rows],
                      [np.nan if r["ca_error_angstrom"] is None else r["ca_error_angstrom"] for r in rows], label=condition)
        confidence_ax.plot([r["target_position"] for r in rows],
                           [np.nan if r["ca_bfactor"] is None else r["ca_bfactor"] for r in rows], label=condition)
    error_ax.set(title=f"{target}: uploaded structures against experimental reference", ylabel="Cα error after rigid fit (Å)")
    confidence_ax.set(ylabel="Cα B value (0–100 expected)", ylim=(0, 100))
    names = ["reference", *conditions, "scoring mask"]
    mask = set(positions)
    coverage = [[int(p in (mask if name == "scoring mask" else mappings[name])) for p in range(1, length + 1)] for name in names]
    coverage_ax.imshow(coverage, aspect="auto", cmap="Greys", vmin=0, vmax=1,
                       extent=(0.5, length + 0.5, len(names) - 0.5, -0.5), interpolation="nearest")
    coverage_ax.set(yticks=range(len(names)), yticklabels=names, xlabel="Target residue position (1-based)")
    for axis in (error_ax, confidence_ax):
        axis.set_xlim(0.5, length + 0.5)
        axis.legend(fontsize=8)
    filename = "per_residue_error_confidence_coverage.png"
    figure.savefig(dest / filename, dpi=180)
    return filename


def audit_structures(root: Path, output_dir: Path, targets: dict) -> dict:
    """Audit all supplied conditions on one frozen residue mask per target.

    Each target specifies ``sequence``, ``reference_path``, ``reference_chain``
    (default A), optional ``reference_mapping``, and ``conditions``. Each
    condition supplies ``model_path``, optional ``model_chain``/``mapping``,
    ``alignment_path`` and ``reported_condition``. Explicit mappings are checked
    for residue identity by the shared metrics module. Ambiguity blocks scoring
    for the whole target, rather than silently removing a condition.
    """
    root, output_dir = Path(root).resolve(), Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    result = {"status": "COMPLETE", "data_kind": "UPLOADED_COORDINATE_AUDIT",
              "rows": [], "artifacts": [], "issues": [], "targets": {}, "inputs": [],
              "metric_definitions": metric_definitions(),
              "comparison_scope": "Static comparison of supplied coordinate files; no causal attribution to MSA strategy.",
              "experimental_controls_status": "UNVERIFIED",
              "unverified_items": ["Complete ColabFold job settings and matched controls", "Rank selection from complete model set",
                                   "pTM from raw score JSON", "Optimized TM-score implementation and full output"],
              "confidence_policy": "Per-residue Cα B values, interpreted as pLDDT only under the uploaded ColabFold model provenance; 0–100 scale. These are not experimental B factors. Means weight each mapped residue equally, not each atom."}
    reported_path = root / "analysis" / "casp15_structure_comparison_metrics.csv"
    reported_rows = []
    if reported_path.is_file():
        with reported_path.open(encoding="utf-8-sig", newline="") as handle:
            reported_rows = list(csv.DictReader(handle))
    for evidence in (reported_path, root / "analysis" / "analysis_summary.md", root / "targets" / "target_selection.md"):
        if evidence.is_file():
            result["inputs"].append({"path": str(evidence), "sha256": _hash(evidence), "role": "uploaded_provenance_document"})
    for target, config in targets.items():
        entry = {"status": "BLOCKED", "inputs": {}, "experimental_controls_status": "UNVERIFIED"}
        result["targets"][target] = entry
        conditions = config.get("conditions", {})
        try:
            if not re.fullmatch(r"[A-Za-z0-9_-]+", target):
                raise ValueError("Target ID must be a safe folder name")
            sequence = config["sequence"]
            if not sequence or set(sequence) - AA:
                raise MappingError("Target sequence must be nonempty canonical protein")
            if not conditions:
                raise MappingError("No uploaded model conditions specified")
            if "reference" in conditions:
                raise ValueError("Condition name 'reference' is reserved")
            ref_path = _path(root, config["reference_path"])
            structures = {"reference": read_ca_structure(ref_path, config.get("reference_chain", "A"), config.get("reference_model_index"))}
            mappings = {"reference": map_to_target(sequence, structures["reference"], config.get("reference_mapping"))}
            mapping_policies = {"reference": "explicit identity-validated" if config.get("reference_mapping") is not None else "unique optimal sequence alignment or exact sequence"}
            entry["reference_header"] = _header(ref_path)
            expected_pdb = CASP15_ASSIGNMENTS.get(target)
            entry["casp_assignment"] = {
                "source": CASP15_TARGET_LIST, "documented_pdb": expected_pdb,
                "uploaded_header_matches_assignment": expected_pdb == entry["reference_header"]["pdb_id"] if expected_pdb else None,
                "chain_selection_basis": "Caller-selected chain, sequence-validated against target; exact mapping saved separately",
            }
            if expected_pdb:
                actual_pdb = entry["reference_header"]["pdb_id"]
                if actual_pdb != expected_pdb:
                    raise MappingError(
                        f"{target}: experimental reference PDB HEADER ID {actual_pdb!r} "
                        f"does not match documented CASP15 assignment {expected_pdb}")
                method = entry["reference_header"]["experimental_method"] or ""
                if "X-RAY DIFFRACTION" not in method.upper():
                    raise MappingError(
                        f"{target}: reference {expected_pdb} lacks the required experimental "
                        "EXPDTA X-RAY DIFFRACTION evidence in its uploaded PDB header")
            for condition, detail in conditions.items():
                path = _path(root, detail["model_path"])
                structures[condition] = read_ca_structure(path, detail.get("model_chain", "A"), detail.get("model_index"))
                mappings[condition] = map_to_target(sequence, structures[condition], detail.get("mapping"))
                mapping_policies[condition] = "explicit identity-validated" if detail.get("mapping") is not None else "unique optimal sequence alignment or exact sequence"
                if detail.get("alignment_path"):
                    alignment = _path(root, detail["alignment_path"])
                    entry["inputs"][condition + "_alignment"] = {"path": str(alignment), "sha256": _hash(alignment)}
            # Freeze intersection before examining any coordinate error or score.
            positions = sorted(set.intersection(*(set(mapping) for mapping in mappings.values())))
            if len(positions) < 3:
                raise MappingError("Fewer than three residues shared across reference and ALL supplied conditions")
            entry.update({"status": "ANALYZED", "target_length": len(sequence), "scored_residue_count": len(positions),
                          "scored_fraction_of_target": len(positions) / len(sequence), "frozen_mask_positions": positions,
                          "coordinate_comparison_status": "COMPARABLE", "mask_policy": "Intersection of sequence-mapped observed Cα positions in reference and every uploaded condition, fixed before fitting; no outlier rejection."})
            for name, structure in structures.items():
                mapped_ids = {row["residue_id"] for row in mappings[name].values()}
                entry["inputs"][name] = {"path": structure["path"], "sha256": structure["sha256"], "chain": structure["chain"],
                    "model_index": structure["model_index"], "coordinate_residue_count": len(structure["residues"]),
                    "mapping_policy": mapping_policies[name], "mapped_count": len(mappings[name]),
                    "unmapped_coordinate_residue_ids": [row["residue_id"] for row in structure["residues"] if row["residue_id"] not in mapped_ids]}
            dest = output_dir / target
            dest.mkdir(parents=True, exist_ok=True)
            mapping_rows = []
            for position, residue in enumerate(sequence, 1):
                row = {"target_position": position, "target_residue": residue, "in_frozen_mask": position in positions}
                for name, mapping in mappings.items():
                    observed = mapping.get(position, {})
                    row.update({f"{name}_{key}": observed.get(key, "") for key in ("chain", "residue_id", "resname", "altloc")})
                mapping_rows.append(row)
            _csv(dest / "residue_mapping.csv", mapping_rows)
            _csv(dest / "frozen_mask.csv", [{"target_position": position} for position in positions])
            result["artifacts"] += [f"{target}/residue_mapping.csv", f"{target}/frozen_mask.csv"]
            reference = np.asarray([mappings["reference"][position]["xyz"] for position in positions])
            profiles = {}
            for condition, detail in conditions.items():
                label = _condition_label(condition, detail)
                reports = [row for row in reported_rows if row.get("Target") == target and row.get("Condition") == label]
                if len(reports) > 1:
                    result["issues"].append(f"{target}/{condition}: multiple uploaded summary rows; comparison withheld")
                uploaded = reports[0] if len(reports) == 1 else {}
                mobile = np.asarray([mappings[condition][position]["xyz"] for position in positions])
                fit, local = kabsch_rmsd(mobile, reference), lddt_ca(mobile, reference)
                errors = dict(zip(positions, fit["distances_angstrom"]))
                rows = [{"target": target, "condition": condition, "target_position": position,
                         "target_residue": residue, "in_frozen_mask": position in errors,
                         "ca_error_angstrom": float(errors[position]) if position in errors else None,
                         "ca_bfactor": mappings[condition].get(position, {}).get("ca_bfactor"),
                         "model_residue_id": mappings[condition].get(position, {}).get("residue_id"),
                         "reference_residue_id": mappings["reference"].get(position, {}).get("residue_id")}
                        for position, residue in enumerate(sequence, 1)]
                profiles[condition] = rows
                bvalues = [row["ca_bfactor"] for row in rows if row["ca_bfactor"] is not None]
                bvalid = all(np.isfinite(value) and 0 <= value <= 100 for value in bvalues)
                bmean = float(np.mean(bvalues)) if bvalues and bvalid else None
                mask_bmean = float(np.mean([mappings[condition][p]["ca_bfactor"] for p in positions])) if bvalid else None
                count = uploaded.get("Mapped_residues")
                same_count = count is None or count == "" or float(count) == len(positions)
                rmsd_check = rounding_comparison(fit["rmsd_angstrom"], uploaded.get("RMSD_CA_A"))
                lddt_check = rounding_comparison(100 * local["score"] if local["score"] is not None else None, uploaded.get("lDDT_CA"))
                mean_check = rounding_comparison(bmean, uploaded.get("pLDDT"))
                metric_row = {"target": target, "condition": condition, "reported_condition": label,
                    "coordinate_comparison_status": "COMPARABLE", "experimental_controls_status": "UNVERIFIED",
                    "target_length": len(sequence), "scored_residue_count": len(positions), "scored_fraction_of_target": len(positions) / len(sequence),
                    "reported_mask_count_matches": same_count if uploaded else None,
                    "rmsd_ca_angstrom": fit["rmsd_angstrom"], "lddt_ca_common_mask": local["score"],
                    "lddt_eligible_reference_pairs": local["eligible_reference_pairs"],
                    "ca_bfactor_residue_mean_all_mapped": bmean, "ca_bfactor_residue_mean_common_mask": mask_bmean,
                    "ca_bfactor_residue_count": len(bvalues), "ca_bfactor_0_100_valid": bvalid,
                    "reported_rmsd_ca_angstrom": rmsd_check["reported"], "rmsd_difference": rmsd_check["difference"],
                    "rmsd_rounding_tolerance": rmsd_check["tolerance"], "rmsd_report_comparison": rmsd_check["status"] if same_count else "NOT_COMPARABLE_MASK_COUNT",
                    "reported_lddt_ca_percent": lddt_check["reported"], "lddt_percent_difference": lddt_check["difference"],
                    "lddt_percent_rounding_tolerance": lddt_check["tolerance"], "lddt_report_comparison": lddt_check["status"] if same_count else "NOT_COMPARABLE_MASK_COUNT",
                    "reported_plddt": mean_check["reported"], "ca_residue_mean_minus_reported_plddt": mean_check["difference"],
                    "plddt_report_comparison": mean_check["status"] if mean_check["status"] != "DIFFERS" else "DIFFERS_AGGREGATION_UNVERIFIED",
                    "reported_ptm": uploaded.get("pTM") or None, "ptm_verification": "UNVERIFIED_NO_RAW_SCORE_JSON",
                    "reported_tm_score": uploaded.get("TM_score") or None, "tm_score_verification": "UNVERIFIED_NO_SCORING_IMPLEMENTATION_OR_FULL_OUTPUT"}
                result["rows"].append(metric_row)
                if not bvalid:
                    result["issues"].append(f"{target}/{condition}: Cα B values outside 0–100; pLDDT interpretation withheld")
                for metric in ("rmsd", "lddt"):
                    state = metric_row[metric + "_report_comparison"]
                    if state in {"DIFFERS", "NOT_COMPARABLE_MASK_COUNT"}:
                        result["issues"].append(f"{target}/{condition}: {metric} summary comparison {state}")
            _csv(dest / "per_residue_profiles.csv", [row for group in profiles.values() for row in group])
            result["artifacts"].append(f"{target}/per_residue_profiles.csv")
            plot = _plots(dest, target, len(sequence), list(conditions), mappings, positions, profiles)
            if plot:
                result["artifacts"].append(f"{target}/{plot}")
            else:
                result["issues"].append(f"{target}: plots unavailable because matplotlib is missing")
        except (MappingError, ValueError, OSError, KeyError, TypeError) as exc:
            entry["status"] = "BLOCKED"
            entry["diagnostic"] = str(exc)
            result["issues"].append(f"{target}: BLOCKED: {exc}")
            result["status"] = "PARTIAL"
    _csv(output_dir / "coordinate_comparison.csv", result["rows"])
    if result["rows"]:
        result["artifacts"].append("coordinate_comparison.csv")
    result["artifacts"].append("result.json")
    (output_dir / "result.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return result
