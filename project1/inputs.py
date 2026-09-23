"""Validated inputs; originals are never edited and upload files contain no placeholders."""
from __future__ import annotations

import copy
import hashlib
import json
import os
import re
from collections import Counter
from pathlib import Path

AA = frozenset("ACDEFGHIKLMNPQRSTVWY")
EXCLUDED = {".git", ".venv", "venv", "env", "tests", "test", "fixtures", "output",
            "outputs", "data", "submissions", "downloads", "__pycache__", ".ipynb_checkpoints"}


def sha256(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    text = json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    if not path.exists() or path.read_text(encoding="utf-8") != text:
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_text(text, encoding="utf-8")
        temporary.replace(path)


def write_text(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not path.exists() or path.read_text(encoding="utf-8") != value:
        path.write_text(value, encoding="utf-8", newline="\n")


def read_fasta(path):
    """All records are candidates; invalid records retain their exact offending sequence."""
    path = Path(path)
    raw = path.read_bytes()
    text = raw.decode("utf-8-sig")
    records, header, lines = [], None, []

    def finish():
        if header is None:
            return
        joined = "".join(lines)
        sequence = "".join(joined.split()).upper()
        ident = header.split()[0] if header.split() else ""
        errors = []
        if not ident:
            errors.append("Empty FASTA identifier")
        if not sequence:
            errors.append("Empty sequence")
        invalid = sorted(set(sequence) - AA)
        if invalid:
            errors.append("Unsupported protein symbols: " + repr("".join(invalid)))
        if ident and not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", ident):
            errors.append("Identifier is not safe for upload filenames; use a documented input override with a safe header")
        records.append({"id": ident, "header": header, "sequence": sequence,
                        "sequence_length": len(sequence), "sequence_sha256": sha256(sequence.encode()),
                        "source_path": str(path.resolve()), "source_sha256": sha256(raw),
                        "normalization": "Derived sequence only: whitespace removed and case uppercased",
                        "normalization_changed": joined != sequence, "errors": errors})

    for line_number, line in enumerate(text.splitlines(), 1):
        if line.startswith(">"):
            finish()
            header, lines = line[1:], []
        elif line.strip():
            if header is None:
                raise ValueError(f"{path}:{line_number}: sequence before first FASTA header")
            lines.append(line)
    finish()
    if not records:
        raise ValueError(f"{path}: no FASTA records")
    return records


def discover_fasta(repo_root, overrides=None, excluded_roots=()):
    root = Path(repo_root).resolve()
    excluded = [Path(p).resolve() for p in excluded_roots]
    paths = []
    for base, dirs, files in os.walk(root):
        dirs[:] = [d for d in dirs if d.lower() not in EXCLUDED and
                   not any((Path(base) / d).resolve().is_relative_to(p) for p in excluded)]
        paths.extend(Path(base) / f for f in sorted(files) if Path(f).suffix.lower() in {".fa", ".fasta", ".faa"})
    # Overrides replace only the named ID, never infer an ID from a filename.
    overrides = overrides or {}
    records, issues = [], []
    for path in sorted(paths):
        try:
            records.extend(r for r in read_fasta(path) if r["id"] not in overrides)
        except (OSError, ValueError, UnicodeError) as exc:
            issues.append(str(exc))
    for ident, value in overrides.items():
        path = Path(value)
        if not path.is_absolute():
            path = root / path
        try:
            candidates = read_fasta(path)
            matches = [r for r in candidates if r["id"] == ident]
            if len(matches) != 1:
                raise ValueError(f"Input override {path} must contain exactly one record with ID {ident}")
            records.extend(matches)
        except (OSError, ValueError, UnicodeError) as exc:
            issues.append(str(exc))
    ids = Counter(r["id"] for r in records)
    seqs = Counter(r["sequence_sha256"] for r in records)
    for r in records:
        if ids[r["id"]] > 1:
            r["errors"].append("Duplicate FASTA identifier")
        if seqs[r["sequence_sha256"]] > 1:
            r["errors"].append("Duplicate target sequence")
        r["validation_status"] = "BLOCKED" if r["errors"] else "VALIDATED"
    return records, issues


def parse_a3m(text, target_sequence):
    """Preserve text for server input; strip lowercase and dots ONLY in analysis rows."""
    headers, rows, pieces = [], [], []
    for line_no, line in enumerate(text.splitlines(), 1):
        if line.startswith(">"):
            if headers:
                rows.append("".join(pieces))
            headers.append(line[1:])
            pieces = []
        elif line.strip():
            if not headers:
                raise ValueError(f"A3M line {line_no}: sequence before header")
            if any(c.isspace() for c in line):
                raise ValueError(f"A3M line {line_no}: whitespace within sequence")
            pieces.append(line)
    if headers:
        rows.append("".join(pieces))
    if not rows or any(not row for row in rows):
        raise ValueError("A3M has no records or an empty record")
    permitted = AA | frozenset(c.lower() for c in AA) | frozenset("Xx-.")
    for i, row in enumerate(rows):
        bad = set(row) - permitted
        if bad:
            raise ValueError(f"A3M row {i + 1}: unsupported symbols {sorted(bad)}")
    match_rows = ["".join(c for c in r if not c.islower() and c != ".") for r in rows]
    # A query containing insertions or gaps is not silently repaired.
    if rows[0] != target_sequence:
        raise ValueError("A3M first sequence must exactly match the submitted target (no gaps or lowercase insertions)")
    if any(len(r) != len(target_sequence) for r in match_rows):
        raise ValueError("A3M inconsistent match-column widths after excluding lowercase insertions and '.'")
    if len(rows) < 2:
        raise ValueError("Query-only A3M is not a real homolog alignment; no DEEPMSA upload generated")
    return {"headers": headers, "rows": rows, "analysis_rows": match_rows,
            "sequence_count_including_query": len(rows), "match_columns": len(target_sequence),
            "sha256": sha256(text.encode("utf-8")), "text": text}


def select_a3m(candidates, explicit=None, service_final=None):
    paths = [Path(p).resolve() for p in candidates]
    if explicit:
        requested = Path(explicit)
        matches = [p for p in paths if str(p) == str(requested.resolve()) or
                   (not requested.is_absolute() and (p.name == str(requested) or p.as_posix().endswith('/' + requested.as_posix())))]
        if len(matches) != 1:
            raise ValueError("Explicit A3M selection must uniquely identify one inventoried candidate")
        return matches[0], "Explicit user configuration selection"
    if service_final:
        if not service_final.get("evidence"):
            raise ValueError("Service final selection requires recorded ranking/final evidence")
        chosen, _ = select_a3m(paths, service_final["path"])
        return chosen, "Service-declared final: " + service_final["evidence"]
    if len(paths) == 1:
        return paths[0], "Only one returned A3M candidate"
    if not paths:
        return None, "WAITING_FOR_DEEPMSA"
    raise ValueError("Multiple A3M candidates: set targets.<ID>.deepmsa.a3m_selection; no automatic ranking assumed")


def server_request(target_id, sequence, seed=42, msa=None):
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed < 2**32:
        raise ValueError("Seed must be a uint32 integer")
    if not sequence or set(sequence) - AA:
        raise ValueError("Invalid protein sequence; no upload generated")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", target_id):
        raise ValueError("Unsafe target ID")
    chain = {"sequence": sequence, "count": 1, "useStructureTemplate": False}
    condition = "AUTO"
    fingerprint = sha256(sequence.encode())[:10]
    if msa is not None:
        parsed = parse_a3m(msa, sequence)
        chain["unpairedMsa"] = msa
        condition = "DEEPMSA"
        fingerprint += "_" + parsed["sha256"][:10]
    request = [{"name": f"p1_{target_id}_{condition}_s{seed}_{fingerprint}", "modelSeeds": [str(seed)],
                "sequences": [{"proteinChain": chain}], "dialect": "alphafoldserver", "version": 1}]
    return request


def verify_pair(auto, deep):
    a, b = copy.deepcopy(auto), copy.deepcopy(deep)
    if "unpairedMsa" in a[0]["sequences"][0]["proteinChain"]:
        raise ValueError("AUTO must omit unpairedMsa")
    msa = b[0]["sequences"][0]["proteinChain"].pop("unpairedMsa")
    parse_a3m(msa, a[0]["sequences"][0]["proteinChain"]["sequence"])
    a[0].pop("name"); b[0].pop("name")
    if a != b:
        raise ValueError("Paired requests differ beyond job identity and MSA input")
    return True
