"""Analyze real, manually uploaded single-chain A3Ms without submitting jobs.

Unlike the submission validator, this analysis accepts a query-only search result:
it documents the lack of homolog evidence instead of inventing a custom MSA.
All statistics use every supplied row. Only the overview image limits its rows.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np

from .metrics import AA, alignment_statistics


def read_uploaded_a3m(path: str | Path, target_sequence: str) -> tuple[list[str], list[str], dict]:
    """Validate A3M query/widths and return insertion-stripped derived rows.

    The optional ColabFold metadata must describe exactly one chain of the
    target length and one copy. Lowercase insertions and dot placeholders are
    removed, never promoted to aligned residues. Originals are never edited.
    """
    if not target_sequence or set(target_sequence) - AA:
        raise ValueError("Target must be a nonempty, ungapped canonical protein sequence")
    names, rows, fragments = [], [], []
    metadata = None
    removed_insertions = 0
    removed_dots = 0

    def finish_record() -> None:
        nonlocal removed_insertions, removed_dots
        if not names:
            return
        raw = "".join(fragments)
        if not raw:
            raise ValueError(f"Empty A3M sequence for row {len(names)}")
        if re.search(r"[^A-Za-z.\-]", raw):
            raise ValueError(f"Unsupported A3M symbols in row {len(names)}")
        removed_insertions += sum("a" <= letter <= "z" for letter in raw)
        removed_dots += raw.count(".")
        row = re.sub(r"[a-z.]", "", raw)
        if set(row) - (AA | frozenset("BXZJUO-")):
            raise ValueError(f"Unsupported A3M amino acids in row {len(names)}")
        if len(row) != len(target_sequence):
            raise ValueError(f"A3M row {len(names)} has width {len(row)}, expected {len(target_sequence)} after removing insertions")
        rows.append(row)
        fragments.clear()

    with Path(path).open(encoding="utf-8-sig") as handle:
        for line_number, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith("#"):
                if names or metadata is not None:
                    raise ValueError("ColabFold metadata is permitted once before the first FASTA record")
                match = re.fullmatch(r"#(\d+)\s+1", line)
                if not match or int(match.group(1)) != len(target_sequence):
                    raise ValueError("ColabFold metadata must match single-chain target length and copy count 1")
                metadata = {"chain_lengths": [int(match.group(1))], "copy_counts": [1]}
            elif line.startswith(">"):
                finish_record()
                name = line[1:].strip()
                if not name:
                    raise ValueError(f"Empty FASTA identifier at line {line_number}")
                names.append(name)
            else:
                if not names:
                    raise ValueError(f"Sequence before FASTA identifier at line {line_number}")
                fragments.append(line)
    finish_record()
    if not rows:
        raise ValueError("A3M contains no FASTA records")
    if rows[0] != target_sequence:
        raise ValueError("A3M query does not exactly match the declared target sequence")
    return names, rows, {
        "colabfold_metadata": metadata,
        "lowercase_insertion_residues_removed": removed_insertions,
        "dot_placeholders_removed": removed_dots,
        "query_coordinate_policy": "Remove lowercase insertion residues and dot placeholders; retain uppercase residues and '-' gaps. Require exact target query and equal final row widths.",
    }


def canonical_conservation_profile(rows: list[str], length: int | None = None) -> list[dict]:
    """Unweighted canonical Shannon entropy and consensus, per alignment column.

    Gaps and noncanonical residues are excluded from the denominator. A column
    with no canonical observations is missing (None), not perfectly conserved.
    Empty rows are useful for the non-query profile of a query-only result.
    """
    width = len(rows[0]) if rows else length
    if width is None or width <= 0 or any(len(row) != width for row in rows):
        raise ValueError("Require equal-width alignment rows or a positive length")
    array = np.frombuffer("".join(rows).encode("ascii"), dtype="S1").reshape(len(rows), width)
    canonical_bytes = [letter.encode("ascii") for letter in sorted(AA)]
    result = []
    for position in range(width):
        counts = np.array([np.count_nonzero(array[:, position] == letter) for letter in canonical_bytes], dtype=int)
        count = int(counts.sum())
        probabilities = counts[counts > 0] / count if count else np.array([])
        entropy = float(-np.sum(probabilities * np.log2(probabilities))) if count else None
        normalized = entropy / math.log2(20) if entropy is not None else None
        result.append({
            "canonical_count": count,
            "shannon_entropy_bits": entropy,
            "shannon_entropy_normalized_log2_20": normalized,
            "conservation_one_minus_normalized_entropy": 1 - normalized if normalized is not None else None,
            "consensus_fraction": float(counts.max() / count) if count else None,
        })
    return result


def _write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _figures(dest: Path, rows: list[str], positions: list[dict], row_metrics: list[dict], status: str, label: str) -> dict:
    from matplotlib.colors import ListedColormap
    from matplotlib.figure import Figure
    from matplotlib.patches import Patch

    length = len(rows[0])
    x = np.arange(1, length + 1)
    fig = Figure(figsize=(12, 7), layout="constrained")
    axes = fig.subplots(2, 2)
    fig.suptitle(f"{label}: {len(rows):,} uploaded rows · {status}", fontsize=12)
    axes[0, 0].plot(x, [p["non_gap_coverage"] for p in positions], label="Non-gap coverage")
    axes[0, 0].plot(x, [p["gap_fraction"] for p in positions], label="Gap fraction")
    axes[0, 0].set(xlabel="Target residue (1-based)", ylabel="Fraction of all rows, query included", ylim=(-0.03, 1.03))
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].plot(x, [p["non_query_conservation_one_minus_normalized_entropy"] for p in positions], label="1 − normalized entropy")
    axes[0, 1].plot(x, [p["non_query_consensus_fraction"] for p in positions], label="Consensus fraction", alpha=0.8)
    axes[0, 1].set(xlabel="Target residue (1-based)", ylabel="Canonical residues in rows after query", ylim=(-0.03, 1.03))
    axes[0, 1].legend(fontsize=8)
    if len(rows) == 1:
        axes[0, 1].text(0.5, 0.5, "Query only: no non-query conservation data", ha="center", transform=axes[0, 1].transAxes, fontsize=9)
    for ax, field, label in ((axes[1, 0], "identity_to_query", "Identity to query (canonical denominator)"),
                              (axes[1, 1], "coverage", "Non-gap coverage / target length")):
        values = [row[field] for row in row_metrics[1:] if row[field] is not None]
        if values:
            ax.hist(values, bins=np.linspace(0, 1, 26), color="#347d96", edgecolor="white")
        else:
            ax.text(0.5, 0.5, "No eligible non-query observations", ha="center", transform=ax.transAxes, fontsize=9)
        ax.set(xlabel=label, ylabel="Rows after query (all eligible rows)", xlim=(0, 1))
    fig.savefig(dest / "alignment_profiles.png", dpi=170)

    # Preserve file order: the query plus evenly spaced row indices. This is a
    # display of at most 80 rows, never a sampled input to any numerical metric.
    chosen = np.unique(np.linspace(0, len(rows) - 1, min(80, len(rows)), dtype=int))
    array = np.frombuffer("".join(rows[int(i)] for i in chosen).encode("ascii"), dtype="S1").reshape(len(chosen), length)
    query = np.frombuffer(rows[0].encode("ascii"), dtype="S1")
    categories = np.full(array.shape, 2, dtype=np.uint8)
    categories[array == b"-"] = 0
    canonical = np.isin(array, [a.encode("ascii") for a in AA])
    categories[(array != b"-") & ~canonical] = 1
    categories[array == query] = 3
    colors = ["#efefef", "#c28a2e", "#9cb6c6", "#1c6f72"]
    fig = Figure(figsize=(12, 5), layout="constrained")
    ax = fig.subplots()
    ax.imshow(categories, cmap=ListedColormap(colors), vmin=-0.5, vmax=3.5, aspect="auto", interpolation="nearest", extent=(0.5, length + 0.5, len(chosen) - 0.5, -0.5))
    tick_positions = np.unique(np.linspace(0, len(chosen) - 1, min(6, len(chosen)), dtype=int))
    ax.set_yticks(tick_positions, [str(int(chosen[i]) + 1) for i in tick_positions])
    ax.set(xlabel="Target residue (1-based)", ylabel="Original file row (query = 1)",
           title=f"{label}: overview of {len(chosen):,}/{len(rows):,} rows\nDeterministic display selection only; statistics use every row")
    ax.legend(handles=[Patch(facecolor=color, label=label) for color, label in zip(colors, ["Gap", "Noncanonical", "Canonical mismatch", "Identical to query"])], loc="upper center", bbox_to_anchor=(0.5, -0.11), ncol=4, fontsize=8)
    fig.savefig(dest / "alignment_overview.png", dpi=170)
    return {"overview_row_indices_zero_based": chosen.tolist(), "overview_max_rows": 80,
            "overview_sampling_policy": "Query plus deterministically spaced file-order rows for visualization only; all numerical summaries use all input rows."}


def analyze_uploaded_msa(path: str | Path, target_sequence: str, output_dir: str | Path) -> dict:
    """Write reproducible local statistics/figures from an uploaded A3M.

    Returns a compact JSON-ready summary; full row/position metrics live in CSV.
    This function never certifies upload readiness or infers the search method
    from sequences. Search provenance remains that of the supplied source file.
    """
    source, dest = Path(path).resolve(), Path(output_dir).resolve()
    names, rows, parse_info = read_uploaded_a3m(source, target_sequence)
    output_names = ["alignment_summary.json", "alignment_positions.csv", "alignment_rows.csv", "query_coordinates.fasta", "alignment_profiles.png", "alignment_overview.png"]
    if source in [dest / name for name in output_names]:
        raise ValueError("Output paths must not overwrite the original input")
    stats = alignment_statistics(rows, neff_max_rows=500)
    positions, row_metrics = stats.pop("positions"), stats.pop("rows")
    all_conservation = canonical_conservation_profile(rows)
    non_query_conservation = canonical_conservation_profile(rows[1:], len(target_sequence))
    for position, all_values, non_query_values in zip(positions, all_conservation, non_query_conservation):
        position.update(all_values)
        position.update({"non_query_" + key: value for key, value in non_query_values.items()})
    distinct_non_query = len(set(rows[1:]) - {target_sequence})
    duplicate_query_rows = sum(row == target_sequence for row in rows[1:])
    status = "QUERY_ONLY" if len(rows) == 1 else "NO_NON_QUERY_DIVERSITY" if distinct_non_query == 0 else "ANALYZED"
    summary = {
        **stats, **parse_info, "status": status,
        "source_path": str(source), "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "analysis_scope": "Descriptive statistics of manually uploaded alignment; not a prediction or a submission-readiness certificate.",
        "distinct_non_query_sequences": distinct_non_query,
        "distinct_non_query_policy": "Number of unique query-coordinate strings among rows after row 0, excluding strings identical to the query; gaps remain part of the string. Not proof of homology.",
        "duplicate_query_rows": duplicate_query_rows,
        "conservation_policy": "Unweighted canonical Shannon entropy H=-sum p*log2(p), normalized by log2(20); conservation=1-H/log2(20). Consensus=max count / canonical count. Exclude gaps and noncanonical amino acids from denominators; no canonical observations => missing. CSV includes all-row and non_query (row 0 excluded) profiles. Exact duplicates retain their input multiplicity.",
        "query_only_policy": "QUERY_ONLY is a valid analysis of an empty search outcome, not homolog evidence or an upload-ready DeepMSA replacement. NO_NON_QUERY_DIVERSITY denotes multiple rows but no distinct string beyond the query.",
        "statistics_row_policy": "All input rows; query included except fields explicitly named non_query or excluding_query. No numerical downsampling.",
        "artifacts": {name: str(dest / name) for name in output_names},
    }
    dest.mkdir(parents=True, exist_ok=True)
    for row, name in zip(row_metrics, names):
        row["sequence_header"] = name
    _write_csv(dest / "alignment_positions.csv", positions)
    _write_csv(dest / "alignment_rows.csv", row_metrics)
    with (dest / "query_coordinates.fasta").open("w", encoding="utf-8", newline="\n") as handle:
        for name, row in zip(names, rows):
            handle.write(f">{name}\n{row}\n")
    try:
        summary.update(_figures(dest, rows, positions, row_metrics, status, Path(path).stem))
        summary["plot_status"] = "AVAILABLE"
    except ImportError:
        summary["plot_status"] = "UNAVAILABLE: install matplotlib to create local figures"
        for name in ("alignment_profiles.png", "alignment_overview.png"):
            summary["artifacts"].pop(name)
    (dest / "alignment_summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return summary
