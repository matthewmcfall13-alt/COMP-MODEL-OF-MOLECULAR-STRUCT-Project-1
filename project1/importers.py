"""Bounded imports of manually downloaded files; never contacts a server.

Source-directory assignment establishes target/condition/run identity. Filenames
associate public-server artifacts within that run, never establish target identity.
All paths returned are absolute. Original bytes are content-addressed and checked
on every rerun. Extraction rejects links, traversal, Windows special paths, archive
bombs, duplicate member paths, and excessive nesting. A blocked import retains its
inventory and previously copied originals for diagnosis.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import math
import os
import re
import shutil
import stat
import tarfile
import tempfile
import zipfile
from collections import defaultdict
from pathlib import Path, PurePosixPath


OUTPUT_DOCUMENTATION = (
    "https://www.ebi.ac.uk/training/online/courses/alphafold/"
    "alphafold-3-and-alphafold-server/alphafold-server-your-gateway-to-alphafold-3/"
    "interpreting-results-from-alphafold-server/"
)
DOCUMENTATION_RETRIEVED = "2026-09-22"
DEFAULT_LIMITS = {
    "max_member_bytes": 512 * 1024 * 1024,
    "max_total_bytes": 2 * 1024 * 1024 * 1024,
    "max_members": 10000,
    "max_depth": 4,
    "max_compression_ratio": 1000,
}


class ImportBlocked(ValueError):
    """Invalid or unsafe input; caller must report a blocked local import."""


def _is_link(path):
    if path.is_symlink():
        return True
    try:
        tag = getattr(path.lstat(), "st_reparse_tag", None)
        return tag in (getattr(stat, "IO_REPARSE_TAG_MOUNT_POINT", -1), getattr(stat, "IO_REPARSE_TAG_SYMLINK", -2))
    except FileNotFoundError:
        return False


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _relative_name(name):
    name = str(name).replace("\\", "/")
    parts = PurePosixPath(name).parts
    reserved = {"CON", "PRN", "AUX", "NUL"} | {
        f"{prefix}{n}" for prefix in ("COM", "LPT") for n in range(1, 10)
    }
    if not parts or name.startswith("/") or any(
        p in ("..", ".") or ":" in p or "\x00" in p
        or p.rstrip(" .") != p or p.split(".")[0].upper() in reserved
        for p in parts
    ):
        raise ImportBlocked(f"Unsafe archive/source path: {name!r}")
    return Path(*parts)


def _safe_destination(root, relative):
    for ancestor in (root, *root.parents):
        if _is_link(ancestor):
            raise ImportBlocked(f"Link/junction in extraction root: {ancestor}")
    root = root.resolve()
    target = root / _relative_name(relative)
    if not target.resolve().is_relative_to(root):
        raise ImportBlocked(f"Path escapes destination: {relative}")
    cursor = target
    while cursor != root:
        if _is_link(cursor):
            raise ImportBlocked(f"Symlink in output path: {cursor}")
        cursor = cursor.parent
    return target


def _archive_kind(path):
    with path.open("rb") as handle:
        head = handle.read(512)
    if head.startswith((b"PK\x03\x04", b"PK\x05\x06", b"PK\x07\x08")):
        return "zip"
    if head.startswith(b"\x1f\x8b"):
        return "gzip"
    if head[257:262] == b"ustar" or path.suffix.lower() == ".tar":
        return "tar"
    return None


def ingest_files(source_dir: Path, destination: Path, overrides=None):
    """Copy and unpack a *manually assigned* directory or file.

    ``overrides`` may override the limits in DEFAULT_LIMITS, not target identity.
    Independent files proceed when another file is rejected. ``files`` includes
    every original and extracted file; archives retain their own inventory entries.
    Size budget counts all original/extracted bytes, including nested containers.
    """
    source_link = _is_link(Path(source_dir))
    source, destination = Path(source_dir).resolve(), Path(destination).resolve()
    limits = dict(DEFAULT_LIMITS)
    if overrides:
        unknown = set(overrides) - set(limits)
        if unknown:
            raise ValueError(f"Unknown extraction limits: {sorted(unknown)}")
        limits.update(overrides)
    if any(not isinstance(value, (int, float)) or value <= 0 for value in limits.values()):
        raise ValueError("Extraction limits must be positive numbers")
    result = {"status": "WAITING_FOR_DOWNLOAD", "source_dir": str(source),
              "destination": str(destination), "originals": [], "files": [],
              "issues": [], "limits": limits}
    if source_link:
        result.update(status="BLOCKED", issues=["Source path is a symlink/junction; use an explicit regular-file/directory path"])
        return result
    if not source.exists():
        return result
    if destination == source or destination.is_relative_to(source):
        raise ValueError("Import destination must be outside the download/source directory")
    destination.mkdir(parents=True, exist_ok=True)
    state = {"bytes": 0, "members": 0}
    visited = set()

    def copy_stream(stream, target, declared_size=None):
        state["members"] += 1
        if state["members"] > limits["max_members"]:
            raise ImportBlocked("Archive/source member count limit exceeded")
        if declared_size is not None and declared_size > limits["max_member_bytes"]:
            raise ImportBlocked(f"Per-member decompression limit exceeded: {target.name}")
        target.parent.mkdir(parents=True, exist_ok=True)
        fd, tmp_name = tempfile.mkstemp(prefix=".import-", dir=target.parent)
        written = 0
        digest = hashlib.sha256()
        try:
            with os.fdopen(fd, "wb") as out:
                while True:
                    data = stream.read(1024 * 1024)
                    if not data:
                        break
                    written += len(data)
                    state["bytes"] += len(data)
                    if written > limits["max_member_bytes"] or state["bytes"] > limits["max_total_bytes"]:
                        raise ImportBlocked("Decompression/copy byte limit exceeded")
                    digest.update(data)
                    out.write(data)
            hexdigest = digest.hexdigest()
            if target.exists():
                if sha256_file(target) != hexdigest:
                    raise ImportBlocked(f"Immutable retained file differs: {target}")
            else:
                # Exclusive creation avoids replacing any pre-existing original.
                with target.open("xb") as out, open(tmp_name, "rb") as inp:
                    shutil.copyfileobj(inp, out)
            return {"path": str(target), "size": written, "sha256": hexdigest}
        finally:
            Path(tmp_name).unlink(missing_ok=True)

    def unpack(entry, depth=0):
        path = Path(entry["path"])
        kind = _archive_kind(path)
        entry["archive_type"] = kind
        if kind is None or entry["sha256"] in visited:
            return
        if depth >= limits["max_depth"]:
            raise ImportBlocked(f"Nested archive depth limit exceeded: {path.name}")
        visited.add(entry["sha256"])
        root = destination / "extracted" / entry["sha256"]
        chain = entry.get("archive_chain", []) + [entry["sha256"]]
        seen_names = set()

        def member(name, stream, size=None, compressed=None):
            clean = _relative_name(name)
            key = clean.as_posix().casefold()
            if key in seen_names:
                raise ImportBlocked(f"Duplicate archive member path: {name}")
            seen_names.add(key)
            if compressed is not None and size and size / max(compressed, 1) > limits["max_compression_ratio"]:
                raise ImportBlocked(f"Compression ratio limit exceeded: {name}")
            child = copy_stream(stream, _safe_destination(root, clean), size)
            child.update({"archive_chain": chain, "member_name": str(name)})
            result["files"].append(child)
            unpack(child, depth + 1)

        if kind == "zip":
            with zipfile.ZipFile(path) as archive:
                infos = archive.infolist()
                if len(infos) + state["members"] > limits["max_members"]:
                    raise ImportBlocked("ZIP member count limit exceeded")
                # Validate every member path and type before extracting anything.
                for info in infos:
                    _relative_name(info.filename)
                    mode = info.external_attr >> 16
                    if stat.S_ISLNK(mode) or (stat.S_IFMT(mode) not in (0, stat.S_IFREG, stat.S_IFDIR)):
                        raise ImportBlocked(f"ZIP link/special member rejected: {info.filename}")
                    if info.flag_bits & 1:
                        raise ImportBlocked("Encrypted ZIP members are unsupported")
                for info in infos:
                    if not info.is_dir():
                        with archive.open(info) as stream:
                            member(info.filename, stream, info.file_size, info.compress_size)
        elif kind == "tar":
            with tarfile.open(path, mode="r:") as archive:
                for info in archive:
                    _relative_name(info.name)
                    if info.isdir():
                        state["members"] += 1
                        if state["members"] > limits["max_members"]:
                            raise ImportBlocked("TAR member count limit exceeded")
                        continue
                    if not info.isfile():
                        raise ImportBlocked(f"TAR link/special member rejected: {info.name}")
                    with archive.extractfile(info) as stream:
                        member(info.name, stream, info.size)
        else:
            name = path.name[:-3] if path.name.lower().endswith(".gz") else path.name + ".uncompressed"
            if path.name.lower().endswith(".tgz"):
                name = path.name[:-4] + ".tar"
            with gzip.open(path, "rb") as stream:
                member(name or "uncompressed", stream)

    sources = []
    if source.is_file():
        sources = [source]
    else:
        for root, dirs, filenames in os.walk(source, followlinks=False):
            for dirname in list(dirs):
                path = Path(root) / dirname
                if _is_link(path):
                    dirs.remove(dirname)
                    result["issues"].append(f"Source directory link rejected: {path}")
            sources.extend(Path(root) / name for name in sorted(filenames))
            if len(sources) > limits["max_members"]:
                result["issues"].append("Source inventory exceeds member count limit; remaining source files were not imported")
                sources = sources[:int(limits["max_members"])]
                break
    for path in sorted(sources):
        if path.name in (".gitkeep", "README.md"):
            continue
        try:
            if _is_link(path) or not path.is_file():
                raise ImportBlocked(f"Source is not a regular file: {path}")
            if path.stat().st_size > limits["max_member_bytes"]:
                raise ImportBlocked(f"Source file exceeds member byte limit: {path.name}")
            digest = sha256_file(path)
            target = _safe_destination(destination / "originals" / digest, path.name)
            with path.open("rb") as stream:
                entry = copy_stream(stream, target, path.stat().st_size)
            if entry["sha256"] != digest:
                raise ImportBlocked(f"Source changed during import: {path}")
            entry.update({"source_relative": path.name if source.is_file() else path.relative_to(source).as_posix(),
                          "archive_chain": []})
            result["originals"].append(dict(entry))
            result["files"].append(entry)
            unpack(entry)
        except (ImportBlocked, OSError, ValueError, EOFError, tarfile.TarError, zipfile.BadZipFile) as exc:
            result["issues"].append(f"{path.name}: {exc}")
    result["status"] = "BLOCKED" if result["issues"] else "IMPORTED" if result["files"] else "WAITING_FOR_DOWNLOAD"
    result["processed_bytes"] = state["bytes"]
    return result


def _read_json(path):
    if Path(path).stat().st_size > 128 * 1024 * 1024:
        raise ValueError("JSON exceeds 128 MiB parser limit; retained without parsing")
    def reject_constant(value):
        raise ValueError(f"Non-finite JSON value: {value}")
    return json.loads(Path(path).read_text(encoding="utf-8-sig"), parse_constant=reject_constant)


def verify_job_request(returned, expected):
    """Compare public server request evidence without substituting defaults."""
    checks = []
    def check(field, want, actual, exists=True):
        status = "UNVERIFIED" if not exists else "VERIFIED" if type(actual) is type(want) and actual == want else "BLOCKED"
        checks.append({"field": field, "status": status, "expected": want, "returned": actual if exists else None})
    if isinstance(expected, list) and len(expected) == 1:
        expected = expected[0]
    if not isinstance(expected, dict):
        raise ValueError("Expected request must describe exactly one public-server job")
    if isinstance(returned, list) and len(returned) == 1:
        returned = returned[0]
    if not isinstance(returned, dict):
        return {"status": "UNVERIFIED", "checks": [{"field": "job_request", "status": "UNVERIFIED", "reason": "Exactly one returned public-server request is required"}]}
    for field in ("dialect", "version", "modelSeeds"):
        check(field, expected.get(field), returned.get(field), field in returned)
    want_entities = expected.get("sequences", [])
    got_entities = returned.get("sequences")
    check("entity_count", len(want_entities), len(got_entities) if isinstance(got_entities, list) else None, isinstance(got_entities, list))
    if isinstance(got_entities, list):
        for index, want_entity in enumerate(want_entities):
            if index >= len(got_entities):
                break
            got_entity = got_entities[index]
            prefix = f"sequences[{index}]"
            if not isinstance(got_entity, dict):
                check(prefix, "entity dictionary", type(got_entity).__name__)
                continue
            check(prefix + ".type", sorted(want_entity), sorted(got_entity))
            for kind, want in want_entity.items():
                got = got_entity.get(kind)
                if not isinstance(got, dict):
                    check(prefix + "." + kind, "entity dictionary", type(got).__name__, kind in got_entity)
                    continue
                for field in ("sequence", "count", "useStructureTemplate"):
                    if field in want:
                        check(prefix + "." + field, want[field], got.get(field), field in got)
                for field in ("templates", "modifications", "glycans", "maxTemplateDate"):
                    if field in want or field in got:
                        check(prefix + "." + field, want.get(field, [] if field != "maxTemplateDate" else None), got.get(field), field in got)
                if "unpairedMsa" in want:
                    # Hash full text: lowercase insertions and row order are controls.
                    hash_text = lambda text: hashlib.sha256(text.encode("utf-8")).hexdigest() if isinstance(text, str) else None
                    check(prefix + ".unpairedMsa_sha256", hash_text(want["unpairedMsa"]), hash_text(got.get("unpairedMsa")), "unpairedMsa" in got)
                else:
                    # In this documented request schema, omission specifies AUTO.
                    if "unpairedMsa" not in got:
                        check(prefix + ".msa_input", "AUTO (unpairedMsa omitted)", "AUTO (unpairedMsa omitted)")
                    else:
                        checks.append({"field": prefix + ".msa_input", "status": "UNVERIFIED",
                                       "reason": "AUTO preparation omitted unpairedMsa, but returned request embeds an alignment; generated-versus-supplied provenance is not established by its presence"})
    statuses = {item["status"] for item in checks}
    return {"status": "BLOCKED" if "BLOCKED" in statuses else "UNVERIFIED" if "UNVERIFIED" in statuses else "VERIFIED", "checks": checks}


def _finite_number(value, minimum=None, maximum=None):
    valid = isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)
    if not valid or (minimum is not None and value < minimum) or (maximum is not None and value > maximum):
        return None
    return float(value)


def _cif_atom_table(path):
    """Keep mmCIF row order so JSON atom arrays can be checked, never assumed."""
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    cif = MMCIF2Dict(str(path))
    def array(key):
        value = cif.get(key, [])
        return value if isinstance(value, list) else [value]
    chains = array("_atom_site.label_asym_id")
    fields = {"label_chain": chains, "auth_chain": array("_atom_site.auth_asym_id"),
              "label_seq_id": array("_atom_site.label_seq_id"), "auth_seq_id": array("_atom_site.auth_seq_id"),
              "insertion_code": array("_atom_site.pdbx_PDB_ins_code"), "atom_name": array("_atom_site.label_atom_id"),
              "residue_name": array("_atom_site.label_comp_id"), "b_iso": array("_atom_site.B_iso_or_equiv")}
    if not chains or any(len(value) != len(chains) for value in fields.values()):
        raise ValueError("Incomplete mmCIF atom mapping fields")
    return [{key: value[i] for key, value in fields.items()} for i in range(len(chains))]


def _verify_model_identity(model_path, expected):
    """Check deposited label-chain/target coordinates, including incomplete models."""
    want = expected[0] if isinstance(expected, list) and len(expected) == 1 else expected
    proteins = [entity["proteinChain"] for entity in want.get("sequences", []) if "proteinChain" in entity]
    if len(proteins) != 1:
        return {"status": "UNVERIFIED", "reason": "Structure identity check supports the declared one-protein study design only"}
    protein = proteins[0]
    try:
        atoms = _cif_atom_table(model_path)
    except (ImportError, OSError, ValueError, KeyError) as exc:
        return {"status": "UNVERIFIED", "reason": f"CIF atom identity could not be parsed: {exc}"}
    codes = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(), "ARNDCQEGHILKMFPSTWYV"))
    residues = {}
    for atom in atoms:
        if atom["label_seq_id"] not in (".", "?"):
            key = (atom["label_chain"], atom["label_seq_id"])
            previous = residues.get(key)
            if previous and previous != atom["residue_name"]:
                return {"status": "BLOCKED", "reason": "Conflicting residue identity at one CIF label position"}
            residues[key] = atom["residue_name"]
    chains = sorted({chain for chain, _ in residues})
    if len(chains) != protein.get("count"):
        return {"status": "BLOCKED", "reason": "CIF polymer chain count differs from prepared protein count", "chains": chains, "expected_count": protein.get("count")}
    missing, mismatches = [], []
    sequence = protein.get("sequence", "")
    for (chain, seq_id), name in residues.items():
        if not seq_id.isdecimal() or not 1 <= int(seq_id) <= len(sequence):
            mismatches.append({"chain": chain, "label_seq_id": seq_id, "reason": "Position outside submitted sequence"})
        elif codes.get(name) != sequence[int(seq_id) - 1]:
            mismatches.append({"chain": chain, "label_seq_id": seq_id, "returned": name, "expected": sequence[int(seq_id) - 1]})
    for chain in chains:
        missing.extend({"chain": chain, "label_seq_id": str(i)} for i in range(1, len(sequence) + 1) if (chain, str(i)) not in residues)
    return {"status": "BLOCKED" if mismatches else "VERIFIED",
            "chains": chains, "observed_residue_count": len(residues), "missing_positions": missing,
            "mismatches": mismatches, "convention": "mmCIF label_seq_id is target position; auth numbering and insertion codes retained separately"}


def _pae_token_mapping(atoms, chains, residue_ids, size):
    """Verify one token per protein residue against CIF label or author identity."""
    if atoms is None or not isinstance(chains, list) or not isinstance(residue_ids, list) or len(chains) != size or len(residue_ids) != size:
        return {"pae_mapping_status": "UNVERIFIED", "pae_mapping_reason": "Token arrays and parsed CIF residue identities are required"}
    tokens = []
    for chain, residue_id in zip(chains, residue_ids):
        valid_id = isinstance(residue_id, (int, str)) and not isinstance(residue_id, bool)
        if not isinstance(chain, str) or not chain or not valid_id or not str(residue_id).isdecimal() or int(residue_id) < 1:
            return {"pae_mapping_status": "UNVERIFIED", "pae_mapping_reason": "Invalid token chain or residue ID"}
        tokens.append((chain, str(int(residue_id))))
    if len(set(tokens)) != size:
        return {"pae_mapping_status": "UNVERIFIED", "pae_mapping_reason": "Duplicate token chain/residue IDs"}
    fields = ("label_chain", "label_seq_id", "auth_chain", "auth_seq_id", "insertion_code", "residue_name")
    residues = {tuple(atom[key] for key in fields): {key: atom[key] for key in fields}
                for atom in atoms if atom["label_seq_id"] not in (".", "?")}
    resolutions = []
    for convention, chain_key, residue_key in (("label", "label_chain", "label_seq_id"), ("author", "auth_chain", "auth_seq_id")):
        mapping = {(row[chain_key], row[residue_key]): key for key, row in residues.items()}
        if len(mapping) != len(residues) or set(tokens) != set(mapping):
            continue
        identities = [mapping[token] for token in tokens]
        resolutions.append((convention, identities))
    if not resolutions:
        return {"pae_mapping_status": "UNVERIFIED", "pae_mapping_reason": "Token chain/residue IDs do not uniquely cover the CIF protein residues"}
    if any(identities != resolutions[0][1] for _, identities in resolutions[1:]):
        return {"pae_mapping_status": "UNVERIFIED", "pae_mapping_reason": "CIF label and author token mappings conflict"}
    return {"pae_mapping_status": "VERIFIED", "pae_mapping_convention": "/".join(convention for convention, _ in resolutions),
            "pae_token_mapping": [dict(residues[key], token_index=index) for index, key in enumerate(resolutions[0][1])]}


def parse_confidence(model_path, summary=None, full=None):
    """Confidence is not accuracy. Residue mean averages within each residue first.

    Atom-to-residue assignment requires CIF row count and label/auth-chain identity
    to agree with the public JSON atom_chain_ids. Without that evidence, the JSON
    atom mean is labelled atom-weighted and no residue mean is manufactured.
    """
    summary = summary if isinstance(summary, dict) else {}
    full = full if isinstance(full, dict) else {}
    out = {"ptm": _finite_number(summary.get("ptm"), 0, 1),
           "ranking_score": _finite_number(summary.get("ranking_score"), -100, 1.5),
           "plddt_units": "0-100", "pae_units": "angstrom",
           "atom_mean_plddt": None, "residue_mean_plddt": None,
           "residue_plddt": [], "atom_mapping": [], "pae": None, "issues": []}
    atoms = None
    try:
        atoms = _cif_atom_table(model_path)
    except (ImportError, OSError, ValueError, KeyError) as exc:
        out["issues"].append(f"CIF atom mapping unavailable: {exc}")
    scores = full.get("atom_plddts")
    chains = full.get("atom_chain_ids")
    if scores is not None:
        valid = isinstance(scores, list) and bool(scores) and all(_finite_number(value, 0, 100) is not None for value in scores)
        if not valid:
            out["issues"].append("Invalid atom_plddts values; no pLDDT aggregate computed")
        else:
            out["atom_mean_plddt"] = sum(scores) / len(scores)
            consistent = atoms is not None and isinstance(chains, list) and len(scores) == len(atoms) == len(chains)
            if consistent:
                consistent = all(str(chain) in (atom["label_chain"], atom["auth_chain"]) for atom, chain in zip(atoms, chains))
            if consistent:
                try:
                    consistent = all(abs(float(atom["b_iso"]) - float(score)) <= 0.11 for atom, score in zip(atoms, scores))
                except (ValueError, TypeError):
                    consistent = False
            if consistent:
                grouped = defaultdict(list)
                for row, score in zip(atoms, scores):
                    mapped = dict(row, plddt=float(score))
                    out["atom_mapping"].append(mapped)
                    # Protein residue identities use the full chain/number/icode key.
                    if row["label_seq_id"] not in (".", "?"):
                        key = tuple(row[key] for key in ("label_chain", "auth_chain", "label_seq_id", "auth_seq_id", "insertion_code", "residue_name"))
                        grouped[key].append(float(score))
                keys = ("label_chain", "auth_chain", "label_seq_id", "auth_seq_id", "insertion_code", "residue_name")
                out["residue_plddt"] = [dict(zip(keys, key), plddt=sum(values) / len(values), atom_count=len(values)) for key, values in grouped.items()]
                if out["residue_plddt"]:
                    out["residue_mean_plddt"] = sum(row["plddt"] for row in out["residue_plddt"]) / len(out["residue_plddt"])
                out["residue_aggregation"] = "Mean of atom pLDDT within each protein residue, then unweighted mean of residue means; mmCIF row-order/chain/B_iso cross-check (0.11 rounding tolerance)"
            else:
                out["issues"].append("Atom-to-residue mapping not verified; residue-weighted pLDDT unavailable")
    pae = full.get("pae")
    token_chains, token_res = full.get("token_chain_ids"), full.get("token_res_ids")
    if pae is not None:
        size = len(pae) if isinstance(pae, list) else 0
        valid = size > 0 and all(isinstance(row, list) and len(row) == size and all(_finite_number(value, 0) is not None for value in row) for row in pae)
        if valid:
            out["pae"] = pae
            out["token_chain_ids"], out["token_res_ids"] = token_chains, token_res
            out.update(_pae_token_mapping(atoms, token_chains, token_res, size))
        else:
            out["issues"].append("Invalid PAE matrix; preserved in original JSON without analysis")
    return out


def import_alphafold(source_dir, destination, expected_request, primary_override=None):
    """Import one explicitly mapped target/condition/run of public-server results.

    ``primary_override`` is {"path": exact retained path or unique basename,
    "rationale": "documented reason"}. It cannot override mismatched job settings.
    Selection occurs solely by server ranking, before any reference scoring.
    """
    result = ingest_files(source_dir, destination)
    result.update({"models": [], "primary_model": None, "primary_selection": None,
                   "returned_requests": [], "returned_alignments": [],
                   "verification": {"status": "UNVERIFIED", "checks": []},
                   "documentation_url": OUTPUT_DOCUMENTATION,
                   "documentation_retrieved": DOCUMENTATION_RETRIEVED})
    if not result["files"]:
        return result
    by_name = defaultdict(list)
    for entry in result["files"]:
        if not entry.get("archive_type"):
            by_name[Path(entry["path"]).name].append(entry)
        if Path(entry["path"]).suffix.lower() in (".a3m", ".a2m", ".sto"):
            result["returned_alignments"].append(entry)
    jsons = {}
    for name, entries in by_name.items():
        if name.lower().endswith(".json"):
            for entry in entries:
                try:
                    jsons[entry["path"]] = _read_json(entry["path"])
                except (OSError, ValueError) as exc:
                    result["issues"].append(f"JSON could not be parsed ({name}): {exc}")
                if name.endswith("_job_request.json"):
                    result["returned_requests"].append(entry)
    requests = result["returned_requests"]
    unique_requests = {item["sha256"]: item for item in requests}
    if len(unique_requests) == 1:
        entry = next(iter(unique_requests.values()))
        returned = jsons.get(entry["path"])
        result["verification"] = verify_job_request(returned, expected_request)
        returned_job = returned[0] if isinstance(returned, list) and len(returned) == 1 else returned
        if isinstance(returned_job, dict) and isinstance(returned_job.get("sequences"), list):
            for index, entity in enumerate(returned_job["sequences"]):
                protein = entity.get("proteinChain", {}) if isinstance(entity, dict) else {}
                msa = protein.get("unpairedMsa") if isinstance(protein, dict) else None
                if isinstance(msa, str) and msa:
                    content = msa.encode("utf-8")
                    digest = hashlib.sha256(content).hexdigest()
                    path = _safe_destination(Path(destination) / "embedded_alignments" / digest, f"returned_entity_{index}.a3m")
                    path.parent.mkdir(parents=True, exist_ok=True)
                    if path.exists() and sha256_file(path) != digest:
                        result["issues"].append(f"Immutable extracted MSA changed: {path}")
                    elif not path.exists():
                        with path.open("xb") as handle:
                            handle.write(content)
                    result["returned_alignments"].append({"path": str(path), "sha256": digest, "size": len(content),
                                                           "source_json": entry["path"], "source_json_sha256": entry["sha256"],
                                                           "field": f"sequences[{index}].proteinChain.unpairedMsa",
                                                           "provenance": "Returned request field; validate query and provenance before analysis"})
    else:
        result["verification"] = {"status": "UNVERIFIED", "checks": [{"field": "job_request", "status": "UNVERIFIED", "reason": "Missing returned job request" if not requests else "Multiple different requests in one run directory; separate downloads by run"}]}
        if len(unique_requests) > 1:
            result["status"] = "BLOCKED"
    model_entries = [entry for entries in by_name.values() for entry in entries if Path(entry["path"]).suffix.lower() in (".cif", ".mmcif")]
    # Identical extracted/plain duplicates represent one retained model identity.
    model_entries = list({(Path(entry["path"]).name, entry["sha256"]): entry for entry in model_entries}.values())
    for entry in sorted(model_entries, key=lambda item: item["path"]):
        model = dict(entry)
        match = re.fullmatch(r"(.+)_model_(\d+)\.(?:cif|mmcif)", Path(entry["path"]).name, re.I)
        model["server_rank"] = int(match.group(2)) if match else None
        summary, full = {}, {}
        if match:
            for suffix, label in (("summary_confidences", "summary"), ("full_data", "full")):
                name = f"{match.group(1)}_{suffix}_{match.group(2)}.json"
                matches = by_name.get(name, [])
                identities = {item["sha256"] for item in matches}
                if len(identities) == 1:
                    model[label + "_path"] = matches[0]["path"]
                    payload = jsons.get(matches[0]["path"], {})
                    if label == "summary":
                        summary = payload
                    else:
                        full = payload
                elif len(identities) > 1:
                    result["issues"].append(f"Ambiguous confidence artifacts: {name}")
        model["confidence"] = parse_confidence(model["path"], summary, full)
        model["ranking_score"] = model["confidence"]["ranking_score"]
        model["identity"] = _verify_model_identity(model["path"], expected_request)
        result["models"].append(model)
    models = result["models"]
    primary, rationale = None, None
    if primary_override is not None:
        if not isinstance(primary_override, dict) or not str(primary_override.get("rationale", "")).strip():
            result["issues"].append("Primary override requires a nonempty rationale")
        else:
            supplied = str(primary_override.get("path", ""))
            matching = [model for model in models if supplied in (model["path"], Path(model["path"]).name)]
            if len(matching) == 1:
                primary, rationale = matching[0], "Explicit recorded resolution: " + primary_override["rationale"]
            else:
                result["issues"].append("Primary override must identify exactly one retained model")
    elif models:
        scores = [model["ranking_score"] for model in models]
        rank_pairs = [(model["server_rank"], model) for model in models]
        unique_ranks = len({rank for rank, _ in rank_pairs}) == len(models) and all(rank is not None for rank, _ in rank_pairs)
        if all(score is not None for score in scores):
            best = max(scores)
            candidates = [model for model in models if model["ranking_score"] == best]
            if len(candidates) == 1:
                primary, rationale = candidates[0], "Maximum public-server ranking_score"
                if unique_ranks and primary["server_rank"] != min(rank for rank, _ in rank_pairs):
                    primary, rationale = None, None
                    result["issues"].append("ranking_score conflicts with documented filename rank; require primary_override rationale")
            elif unique_ranks:
                primary = min(candidates, key=lambda model: model["server_rank"])
                rationale = "Tied ranking_score resolved by documented public-server filename rank"
        elif unique_ranks:
            primary = min(models, key=lambda model: model["server_rank"])
            rationale = "Documented public-server *_model_N.cif rank (lower is better)"
        if primary is None:
            result["issues"].append("Primary model ranking missing or ambiguous; set primary_override with rationale")
    if primary and primary_override is None and primary.get("server_rank") not in (None, 0):
        result["issues"].append("Highest-ranked server model (rank 0) is absent; download complete results or record a primary_override rationale")
        primary, rationale = None, None
    if primary:
        result["primary_model"], result["primary_selection"] = primary["path"], rationale
    for model in models:
        result["verification"]["checks"].append({"field": "structure_identity", "path": model["path"], **model["identity"]})
    identity_states = {model["identity"]["status"] for model in models}
    if "BLOCKED" in identity_states:
        result["verification"]["status"] = "BLOCKED"
    elif "UNVERIFIED" in identity_states and result["verification"]["status"] == "VERIFIED":
        result["verification"]["status"] = "UNVERIFIED"
    if not models:
        result["issues"].append("No returned CIF models found; original download retained")
    result["automatic_msa_status"] = "AVAILABLE_CANDIDATES_REQUIRE_VALIDATION" if result["returned_alignments"] else "UNAVAILABLE: no returned alignment artifact; structural analysis can proceed independently"
    if result["verification"]["status"] == "BLOCKED" or result["status"] == "BLOCKED":
        result["status"] = "BLOCKED"
    elif result["issues"] or not primary:
        result["status"] = "IMPORTED_NEEDS_REVIEW"
    elif result["verification"]["status"] == "VERIFIED":
        result["status"] = "VALIDATED"
    else:
        result["status"] = "IMPORTED_UNVERIFIED"
    result["comparison_eligible"] = result["status"] == "VALIDATED"
    return result
