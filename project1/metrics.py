"""Local alignment and coordinate analysis; never invokes a prediction service.

Accuracy is distinct from AlphaFold confidence. Missing scores are None, not zero.
Reference contacts and a single paired residue mask are defined before scoring.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
from typing import Any

import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB import MMCIFParser, PDBParser, is_aa
from Bio.SeqUtils import seq1, seq3


AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
TMSCORE_SOURCE = "https://github.com/pylelab/USalign/blob/master/TMscore.cpp"
LDDT_SOURCE = "https://openstructure.org/docs/2.11/mol/alg/lddt/"


class MappingError(ValueError):
    """An unambiguous, sequence-consistent coordinate map is required."""


def _write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def alignment_statistics(
    rows: list[str], output_dir: str | Path | None = None,
    neff_max_rows: int = 500, neff_identity: float = 0.8,
    neff_min_coverage: float = 0.5,
) -> dict:
    """Analyze an A3M-derived *query-coordinate* copy, query included as row 0.

    Lowercase insertions and dots must already have been removed by an A3M-aware
    parser. '-' is the only gap. Ambiguous amino acids are counted in coverage,
    but excluded from BLOSUM62 and identity denominators. The query is canonical.
    """
    if not rows or not rows[0] or any(len(row) != len(rows[0]) for row in rows):
        raise ValueError("Alignment must contain equal-width, nonempty query-coordinate rows")
    if set(rows[0]) - AA:
        raise ValueError("Query row must contain only 20 canonical amino acids without gaps")
    allowed = AA | frozenset("BXZJUO-")
    for index, row in enumerate(rows):
        if set(row) - allowed:
            raise ValueError(f"Row {index + 1}: lowercase insertions or unsupported symbols remain")
    if not 0 <= neff_identity <= 1 or not 0 <= neff_min_coverage <= 1:
        raise ValueError("Neff identity and coverage thresholds must lie in [0, 1]")
    arr = np.array([list(row) for row in rows], dtype="U1")
    n_rows, length = arr.shape
    canonical = np.isin(arr, list(AA))
    gaps = arr == "-"
    matrix = substitution_matrices.load("BLOSUM62")
    total_score = 0.0
    total_pairs = 0
    position_rows = []
    for index in range(length):
        letters, counts = np.unique(arr[canonical[:, index], index], return_counts=True)
        column_pairs = 0
        column_score = 0.0
        for ai, (a, count_a) in enumerate(zip(letters, counts)):
            same_pairs = int(count_a) * (int(count_a) - 1) // 2
            column_score += same_pairs * float(matrix[a, a])
            column_pairs += same_pairs
            for b, count_b in zip(letters[ai + 1:], counts[ai + 1:]):
                pairs = int(count_a) * int(count_b)
                column_score += pairs * float(matrix[a, b])
                column_pairs += pairs
        total_score += column_score
        total_pairs += column_pairs
        position_rows.append({
            "target_position": index + 1, "query_residue": rows[0][index],
            "non_gap_coverage": float(np.mean(~gaps[:, index])),
            "canonical_coverage": float(np.mean(canonical[:, index])),
            "gap_fraction": float(np.mean(gaps[:, index])),
            "eligible_unordered_pairs": column_pairs,
            "blosum62_sum": column_score,
        })
    per_row = []
    for index, row in enumerate(rows):
        denominator = int(canonical[index].sum())
        matches = int(((arr[index] == arr[0]) & canonical[index]).sum())
        per_row.append({
            "row_index": index, "is_query": index == 0,
            "coverage": float(np.mean(~gaps[index])),
            "identity_to_query": matches / denominator if denominator else None,
            "identity_denominator": denominator,
        })
    unique = len(set(rows))
    neff = None
    neff_status = "UNAVAILABLE_EXACT_ROW_LIMIT"
    if n_rows <= neff_max_rows:
        neighbors = np.zeros(n_rows, dtype=int)
        # Bounded O(N^2 L); never silently subsample a large MSA.
        for i in range(n_rows):
            for j in range(i, n_rows):
                shared = canonical[i] & canonical[j]
                overlap = int(shared.sum())
                identity = float(np.mean(arr[i, shared] == arr[j, shared])) if overlap else 0.0
                if overlap / length >= neff_min_coverage and identity >= neff_identity:
                    neighbors[i] += 1
                    if j != i:
                        neighbors[j] += 1
        # Rows failing the coverage threshold have weight 0, even for self.
        neff = float(sum(1.0 / count for count in neighbors if count))
        neff_status = "EXACT"
    identities = [row["identity_to_query"] for row in per_row[1:] if row["identity_to_query"] is not None]
    summary = {
        "query_length": length, "sequence_count_including_query": n_rows,
        "homolog_row_count": n_rows - 1, "unique_query_coordinate_rows": unique,
        "duplicate_rows_beyond_first": n_rows - unique,
        "duplicate_fraction": (n_rows - unique) / n_rows,
        "gap_fraction": float(gaps.mean()),
        "mean_identity_to_query_excluding_query": float(np.mean(identities)) if identities else None,
        "blosum62_sum_of_pairs": total_score,
        "blosum62_eligible_pair_denominator": total_pairs,
        "blosum62_normalized_sum_of_pairs": total_score / total_pairs if total_pairs else None,
        "blosum62_policy": "Unordered pairs across columns, query included; gaps and noncanonical residues excluded. Matrix-score units per eligible residue pair; not a 0-1 score.",
        "identity_policy": "Identical canonical residues / canonical residues in each row; query excluded from reported identity mean.",
        "duplicate_policy": "Exact duplicate query-coordinate strings, including gap positions; insertions not considered.",
        "neff": neff, "neff_status": neff_status,
        "neff_definition": f"Query included; pair identity >= {neff_identity} on jointly canonical columns AND joint coverage/query_length >= {neff_min_coverage}; sum_i 1/neighbors_i; rows with no eligible neighbor weight 0.",
        "neff_max_rows": neff_max_rows, "neff_approximation": False,
        "search_evalues": None, "search_pvalues": None,
        "positions": position_rows, "rows": per_row,
    }
    if output_dir is not None:
        dest = Path(output_dir)
        dest.mkdir(parents=True, exist_ok=True)
        _write_json(dest / "alignment_statistics.json", summary)
        _write_csv(dest / "alignment_positions.csv", position_rows, list(position_rows[0]))
        _write_csv(dest / "alignment_rows.csv", per_row, list(per_row[0]))
        try:
            from matplotlib.figure import Figure
            fig = Figure(figsize=(9, 3.5), layout="constrained")
            ax = fig.subplots()
            ax.plot(range(1, length + 1), [x["non_gap_coverage"] for x in position_rows], label="Non-gap coverage")
            ax.plot(range(1, length + 1), [x["gap_fraction"] for x in position_rows], label="Gap fraction")
            ax.set(xlabel="Target residue position (1-based)", ylabel="Fraction of rows (query included)", ylim=(-0.02, 1.02))
            ax.legend()
            fig.savefig(dest / "alignment_coverage.png", dpi=180)
        except ImportError:
            summary["plot_status"] = "UNAVAILABLE: install matplotlib"
    return summary


def read_ca_structure(path: str | Path, chain: str, model_index: int | None = None) -> dict:
    """Read protein Cα atoms; preserve author chain/number/insertion identifiers."""
    path = Path(path)
    parser = MMCIFParser(QUIET=True) if path.suffix.lower() in {".cif", ".mmcif"} else PDBParser(QUIET=True)
    structure = parser.get_structure("coordinates", str(path))
    models = list(structure.get_models())
    if len(models) != 1 and model_index is None:
        raise MappingError(f"{path.name}: {len(models)} coordinate models; select model_index explicitly")
    index = 0 if model_index is None else model_index
    if index < 0 or index >= len(models):
        raise MappingError(f"Invalid model index {index} for {path.name}")
    if chain not in models[index]:
        raise MappingError(f"{path.name}: chain {chain!r} absent; present: {[c.id for c in models[index]]}")
    residues = []
    missing_ca = []
    for residue in models[index][chain]:
        if not is_aa(residue, standard=False):
            continue
        number, insertion = int(residue.id[1]), residue.id[2].strip()
        identifier = f"{number}{insertion}"
        amino_acid = seq1(residue.resname, custom_map={"MSE": "M", "SEC": "U", "PYL": "O"})
        entry = {"chain": chain, "number": number, "insertion_code": insertion,
                 "residue_id": identifier, "resname": residue.resname, "amino_acid": amino_acid}
        if "CA" not in residue:
            missing_ca.append(entry)
            continue
        atom = residue["CA"]
        entry.update({"xyz": np.asarray(atom.coord, dtype=float), "ca_bfactor": float(atom.bfactor),
                      "altloc": atom.get_altloc().strip(), "occupancy": atom.get_occupancy()})
        residues.append(entry)
    if not residues:
        raise MappingError(f"{path.name}: selected chain has no protein C-alpha coordinates")
    ids = [row["residue_id"] for row in residues]
    if len(ids) != len(set(ids)):
        raise MappingError("Duplicate author residue identifiers require an explicit cleaned coordinate copy")
    return {"path": str(path.resolve()), "sha256": _sha(path), "chain": chain,
            "model_index": index, "residues": residues, "missing_ca": missing_ca,
            "sequence": "".join(row["amino_acid"] for row in residues)}


def map_to_target(target: str, structure: dict, explicit_mapping: dict | None = None) -> dict[int, dict]:
    """Return target-position -> observed residue, blocking ambiguous alignments.

    Explicit mapping values are author residue IDs, e.g. {1: '42', 2: '42A'}.
    Omitted target positions remain unscored. Identity is still checked.
    """
    residues = structure["residues"]
    if explicit_mapping is not None:
        lookup = {r["residue_id"]: r for r in residues}
        result = {}
        used = set()
        for raw_position, raw_identifier in explicit_mapping.items():
            position, identifier = int(raw_position), str(raw_identifier)
            if not 1 <= position <= len(target) or identifier not in lookup or identifier in used:
                raise MappingError(f"Invalid/repeated explicit mapping: target {position} -> residue {identifier}")
            residue = lookup[identifier]
            if residue["amino_acid"] != target[position - 1]:
                raise MappingError(f"Explicit mapping sequence mismatch at target position {position}")
            result[position] = residue
            used.add(identifier)
        if not result:
            raise MappingError("Explicit residue mapping is empty")
        return result
    observed = structure["sequence"]
    if observed == target:
        return {i + 1: r for i, r in enumerate(residues)}
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -3.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -0.5
    alignments = iter(aligner.align(target, observed))
    first = next(alignments)
    if next(alignments, None) is not None:
        raise MappingError("Ambiguous optimal sequence-to-coordinate mapping; provide explicit mapping")
    result = {}
    for (t0, t1), (r0, r1) in zip(first.aligned[0], first.aligned[1]):
        for ti, ri in zip(range(int(t0), int(t1)), range(int(r0), int(r1))):
            if target[ti] != observed[ri]:
                raise MappingError(f"Sequence mismatch at target {ti + 1}; construct differences require explicit mapping")
            result[ti + 1] = residues[ri]
    if not result:
        raise MappingError("No sequence-identical residue correspondences")
    return result


def kabsch_rmsd(mobile: np.ndarray, reference: np.ndarray) -> dict:
    """Least-squares rigid fit without reflection, outlier rejection, or pruning."""
    mobile, reference = np.asarray(mobile, dtype=float), np.asarray(reference, dtype=float)
    if mobile.shape != reference.shape or mobile.ndim != 2 or mobile.shape[1] != 3 or len(mobile) < 3:
        raise ValueError("At least three paired 3D points are required")
    if not np.isfinite(mobile).all() or not np.isfinite(reference).all():
        raise ValueError("Coordinates must be finite")
    m_center, r_center = mobile.mean(axis=0), reference.mean(axis=0)
    u, _, vt = np.linalg.svd((mobile - m_center).T @ (reference - r_center))
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    translation = r_center - m_center @ rotation
    fitted = mobile @ rotation + translation
    distances = np.linalg.norm(fitted - reference, axis=1)
    return {"rmsd_angstrom": float(np.sqrt(np.mean(distances ** 2))),
            "rotation_row_vector": rotation, "translation": translation,
            "fitted": fitted, "distances_angstrom": distances}


def lddt_ca(mobile: np.ndarray, reference: np.ndarray, radius: float = 15.0) -> dict:
    """Cα distance-difference test on the supplied fixed mask, not all-atom lDDT.

    Distinct-residue unordered reference pairs at distance <15 Å are eligible;
    conserved distances have absolute error strictly <0.5,1,2,4 Å. No fitting.
    A NaN mobile coordinate counts as an unconserved reference contact.
    """
    mobile, reference = np.asarray(mobile, float), np.asarray(reference, float)
    if mobile.shape != reference.shape or mobile.ndim != 2 or mobile.shape[1] != 3:
        raise ValueError("Paired Nx3 arrays are required")
    if not np.isfinite(reference).all() or radius <= 0:
        raise ValueError("Reference coordinates must be finite and radius positive")
    ref_dist = np.linalg.norm(reference[:, None] - reference[None, :], axis=2)
    mob_dist = np.linalg.norm(mobile[:, None] - mobile[None, :], axis=2)
    eligible = np.triu(ref_dist < radius, k=1)
    errors = np.abs(mob_dist[eligible] - ref_dist[eligible])
    conserved = [int(np.count_nonzero(errors < threshold)) for threshold in (0.5, 1.0, 2.0, 4.0)]
    pair_count = len(errors)
    return {"score": sum(conserved) / (4 * pair_count) if pair_count else None,
            "eligible_reference_pairs": pair_count, "conserved_by_threshold": conserved,
            "thresholds_angstrom": [0.5, 1.0, 2.0, 4.0], "radius_angstrom": radius,
            "scope": "lDDT-Calpha on supplied common residue mask; no stereochemical checks",
            "source": LDDT_SOURCE}


def parse_tmscore_output(stdout: str, reference_length: int, mask_count: int) -> dict:
    """Parse official TMscore full output; preserve its explicit normalization."""
    if reference_length < mask_count or mask_count < 3:
        raise ValueError("Reference normalization length must be >= common mask size >= 3")
    tm_lines = re.findall(r"^TM-score\s*=\s*([0-9.eE+-]+)\s*\(([^\n]+)\)", stdout, re.MULTILINE)
    tm_score = None
    for score, explanation in tm_lines:
        if "user-specified" in explanation:
            length = re.search(r"LN?\s*=\s*([0-9.]+)", explanation)
            if not length or not math.isclose(float(length.group(1)), reference_length):
                raise ValueError("TMscore output normalization does not match declared reference length")
            tm_score = float(score)
    if tm_score is None:
        raise ValueError("TMscore output lacks user-specified reference-length score (-l)")
    gdt_match = re.search(r"^GDT-TS-score\s*=\s*([0-9.eE+-]+)", stdout, re.MULTILINE)
    native_len = re.search(r"^Structure2:.*?Length=\s*(\d+)", stdout, re.MULTILINE)
    common = re.search(r"Number of residues in common=\s*(\d+)", stdout)
    if not native_len or int(native_len.group(1)) != mask_count or not common or int(common.group(1)) != mask_count:
        raise ValueError("TMscore did not preserve the entire frozen residue correspondence")
    gdt_native = float(gdt_match.group(1)) if gdt_match else None
    if not 0 <= tm_score <= 1 or (gdt_native is not None and not 0 <= gdt_native <= 1):
        raise ValueError("TMscore returned a score outside [0,1]")
    return {"tm_score": tm_score, "gdt_ts": gdt_native * mask_count / reference_length if gdt_native is not None else None,
            "gdt_ts_native_mask_normalized": gdt_native,
            "reference_normalization_length": reference_length,
            "gdt_normalization": "Official optimized GDT counts / declared reference length (rescaled from native mask length)",
            "status": "AVAILABLE", "source": TMSCORE_SOURCE}


def run_tmscore(model_path: Path, reference_path: Path, reference_length: int,
                mask_count: int, executable: str | Path | None = None) -> dict:
    """Optional CPU scoring wrapper. Never installs or downloads an executable."""
    requested = str(executable) if executable else "TMscore"
    resolved = shutil.which(requested)
    if not resolved and Path(requested).is_file():
        resolved = str(Path(requested).resolve())
    requirement = "Provide an official pylelab/USalign TMscore executable (not TMalign) through tmscore_executable; see " + TMSCORE_SOURCE
    if not resolved:
        return {"status": "UNAVAILABLE", "tm_score": None, "gdt_ts": None, "requirement": requirement}
    command = [resolved, str(model_path), str(reference_path), "-l", str(reference_length)]
    try:
        completed = subprocess.run(command, capture_output=True, text=True, timeout=120, check=True)
        parsed = parse_tmscore_output(completed.stdout, reference_length, mask_count)
        parsed.update({"command": command, "executable_sha256": _sha(Path(resolved)),
                       "stdout": completed.stdout, "stderr": completed.stderr})
        return parsed
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        return {"status": "UNAVAILABLE", "tm_score": None, "gdt_ts": None,
                "diagnostic": str(exc), "requirement": requirement}


def _write_ca_pdb(path: Path, positions: list[int], sequence: str, xyz: np.ndarray) -> None:
    lines = []
    for serial, (position, coords) in enumerate(zip(positions, xyz), start=1):
        residue = seq3(sequence[position - 1]).upper()
        x, y, z = coords
        lines.append(f"ATOM  {serial:5d}  CA  {residue:3s} A{position:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}           C  ")
    path.write_text("\n".join(lines + ["TER", "END", ""]), encoding="ascii")


def analyze_pair(
    target_sequence: str, reference_path: str | Path, reference_chain: str,
    auto_path: str | Path, deep_path: str | Path, output_dir: str | Path,
    auto_chain: str = "A", deep_chain: str = "A",
    reference_mapping: dict | None = None, auto_mapping: dict | None = None,
    deep_mapping: dict | None = None, reference_length: int | None = None,
    tmscore_executable: str | Path | None = None,
    reference_model_index: int | None = None,
) -> dict:
    """Score preselected paired primary structures against a confirmed reference.

    Caller must establish reference provenance, target identity and matched job
    controls. This function neither chooses primary models nor confirms a PDB
    assignment. TM/GDT require an explicitly confirmed reference_length; target
    length is not silently substituted for the full experimental-chain length.
    """
    dest = Path(output_dir)
    dest.mkdir(parents=True, exist_ok=True)
    normalization = None if reference_length is None else int(reference_length)
    result = {"status": "BLOCKED", "target_length": len(target_sequence),
              "reference_normalization_length": normalization,
              "normalization_policy": "UNRESOLVED: provide confirmed full experimental-chain reference_length" if reference_length is None else "Explicit declared experimental reference length",
              "reference_assignment": "Must be confirmed by caller; sequence matching alone is not confirmation",
              "metric_definitions": metric_definitions(), "metrics": {}}
    try:
        if not target_sequence or set(target_sequence) - AA:
            raise MappingError("Target sequence must be nonempty canonical protein")
        structures = {
            "reference": read_ca_structure(reference_path, reference_chain, reference_model_index),
            "AUTO": read_ca_structure(auto_path, auto_chain),
            "DEEPMSA": read_ca_structure(deep_path, deep_chain),
        }
        mappings = {
            "reference": map_to_target(target_sequence, structures["reference"], reference_mapping),
            "AUTO": map_to_target(target_sequence, structures["AUTO"], auto_mapping),
            "DEEPMSA": map_to_target(target_sequence, structures["DEEPMSA"], deep_mapping),
        }
        positions = sorted(set(mappings["reference"]) & set(mappings["AUTO"]) & set(mappings["DEEPMSA"]))
        if len(positions) < 3 or (normalization is not None and normalization < len(positions)):
            raise MappingError("Require >=3 common C-alpha residues and normalization length >= scored count")
        mapping_rows = []
        for position in range(1, len(target_sequence) + 1):
            row = {"target_position": position, "target_residue": target_sequence[position - 1], "in_frozen_mask": position in positions}
            for name, mapping in mappings.items():
                residue = mapping.get(position, {})
                row.update({f"{name}_{field}": residue.get(field, "") for field in ("chain", "number", "insertion_code", "residue_id", "resname", "altloc")})
            mapping_rows.append(row)
        _write_csv(dest / "residue_mapping.csv", mapping_rows, list(mapping_rows[0]))
        mask_rows = [{"target_position": position} for position in positions]
        _write_csv(dest / "frozen_mask.csv", mask_rows, ["target_position"])
        coords = {name: np.array([mapping[position]["xyz"] for position in positions]) for name, mapping in mappings.items()}
        for name, xyz in coords.items():
            _write_ca_pdb(dest / f"{name}_common_mask.pdb", positions, target_sequence, xyz)
        provenance = {}
        for name, structure in structures.items():
            provenance[name] = {key: value for key, value in structure.items() if key != "residues"}
            mapped_ids = {r["residue_id"] for r in mappings[name].values()}
            provenance[name]["unmapped_coordinate_residues"] = [r["residue_id"] for r in structure["residues"] if r["residue_id"] not in mapped_ids]
            provenance[name]["mapping_policy"] = "explicit" if {"reference": reference_mapping, "AUTO": auto_mapping, "DEEPMSA": deep_mapping}[name] is not None else "unique optimal sequence alignment (or exact sequence)"
        result.update({"status": "ANALYZED", "inputs": provenance, "scored_residue_count": len(positions),
                       "scored_fraction_of_target": len(positions) / len(target_sequence),
                       "scored_fraction_of_reference_normalization": len(positions) / normalization if normalization is not None else None,
                       "frozen_mask": "frozen_mask.csv", "mapping": "residue_mapping.csv"})
        fit_errors = {}
        for condition in ("AUTO", "DEEPMSA"):
            fit = kabsch_rmsd(coords[condition], coords["reference"])
            local = lddt_ca(coords[condition], coords["reference"])
            if normalization is None:
                standard = {"status": "UNAVAILABLE", "tm_score": None, "gdt_ts": None,
                            "requirement": "Set reference.reference_length to the confirmed full experimental-chain normalization length; target length is not assumed to equal reference length"}
            else:
                standard = run_tmscore(dest / f"{condition}_common_mask.pdb", dest / "reference_common_mask.pdb", normalization, len(positions), tmscore_executable)
            if "stdout" in standard:
                (dest / f"{condition}_TMscore.txt").write_text(standard.pop("stdout"), encoding="utf-8")
            fit_errors[condition] = fit["distances_angstrom"]
            result["metrics"][condition] = {"rmsd_angstrom": fit["rmsd_angstrom"],
                "lddt_ca_common_mask": local["score"], "lddt_details": local,
                "tm_score": standard["tm_score"], "gdt_ts": standard["gdt_ts"], "standard_tool": standard,
                "molprobity": None, "molprobity_status": "UNAVAILABLE: requires separately supported MolProbity/Phenix geometry analysis",
                "rigid_transform": {"rotation_row_vector": fit["rotation_row_vector"].tolist(), "translation": fit["translation"].tolist()}}
        error_rows = [{"target_position": position, "AUTO_ca_error_angstrom": float(fit_errors["AUTO"][i]),
                       "DEEPMSA_ca_error_angstrom": float(fit_errors["DEEPMSA"][i])} for i, position in enumerate(positions)]
        _write_csv(dest / "ca_position_errors.csv", error_rows, list(error_rows[0]))
        result["paired_differences_DEEPMSA_minus_AUTO"] = {
            metric: result["metrics"]["DEEPMSA"][metric] - result["metrics"]["AUTO"][metric]
            if result["metrics"]["DEEPMSA"][metric] is not None and result["metrics"]["AUTO"][metric] is not None else None
            for metric in ("rmsd_angstrom", "lddt_ca_common_mask", "tm_score", "gdt_ts")}
        result["metric_coverage"] = {metric: ("AVAILABLE" if result["metrics"]["AUTO"][metric] is not None and result["metrics"]["DEEPMSA"][metric] is not None else "UNAVAILABLE")
                                     for metric in ("rmsd_angstrom", "lddt_ca_common_mask", "tm_score", "gdt_ts", "molprobity")}
        # Coordinate copies carry the saved target-position mapping. No sequence
        # re-alignment or distance-based pruning is permitted in the overlay.
        script = ["# Common-mask CA-only primary-model/reference overlay; see residue_mapping.csv.",
                  "# Least-squares fit uses every frozen pair; no cutoffDistance pruning."]
        for name in ("reference", "AUTO", "DEEPMSA"):
            script.append(f'open "{(dest / (name + "_common_mask.pdb")).resolve().as_posix()}"')
        script += ["align #2@CA toAtoms #1@CA matchNumbering true matchAtomNames true",
                   "align #3@CA toAtoms #1@CA matchNumbering true matchAtomNames true",
                   "color #1 gray", "color #2 blue", "color #3 orange", "show atoms", "view"]
        (dest / "overlay.cxc").write_text("\n".join(script) + "\n", encoding="utf-8")
        try:
            from matplotlib.figure import Figure
            fig = Figure(figsize=(9, 3.5), layout="constrained")
            ax = fig.subplots()
            for condition, distances in fit_errors.items():
                ax.plot(positions, distances, label=condition)
            ax.set(xlabel="Target residue position (frozen common mask)", ylabel="Cα error after rigid fit (Å)")
            ax.legend()
            fig.savefig(dest / "paired_ca_errors.png", dpi=180)
        except ImportError:
            result["plot_status"] = "UNAVAILABLE: install matplotlib"
    except (MappingError, ValueError, OSError, KeyError) as exc:
        result["status"] = "BLOCKED"
        result["metrics"] = {}
        result.pop("paired_differences_DEEPMSA_minus_AUTO", None)
        result.pop("metric_coverage", None)
        result["diagnostic"] = str(exc)
    _write_json(dest / "metric_summary.json", result)
    return result


def metric_definitions() -> dict:
    return {
        "rmsd_angstrom": "sqrt(mean squared C-alpha distance) after Kabsch rigid fit, determinant +1, same paired mask, no pruning; Angstrom.",
        "lddt_ca_common_mask": "Unordered reference C-alpha pairs <15 Angstrom; mean fractions with absolute distance error strictly <0.5,1,2,4 Angstrom; supplied frozen common mask only; 0-1. Missing reference positions outside mask excluded and coverage reported.",
        "tm_score": "Official optimized TMscore, fixed author-index correspondence in derived renumbered copies, -l declared full reference normalization length. Missing positions contribute no numerator and stay in normalization. UNAVAILABLE without validated executable.",
        "gdt_ts": "Official TMscore optimized GDT-TS count average at 1,2,4,8 Angstrom divided by declared reference length; rescale program mask-length-normalized output. 0-1, not percent and not one-fit threshold fraction.",
        "molprobity": "Separate all-atom geometry assessment; unavailable without supported external tool; never substitute accuracy/confidence scores.",
        "sources": [TMSCORE_SOURCE, LDDT_SOURCE, "presentation_work/source_audit_extracts/2_Sequence_alignment_v7.pdf.txt pp.102-103 (parent workspace)",
                    "https://www.cgl.ucsf.edu/chimerax/docs/user/commands/align.html"],
    }
