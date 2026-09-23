"""Evidence-limited review packages; originals are only ever read.

Tests are synthetic software checks, never study observations. Explicit optional
``status['review_answers']`` entries must identify real evidence files. A prepared
request, an existing output directory, and a passing test do not establish a
server submission or a scientific result.
"""
from __future__ import annotations

import csv
from datetime import datetime, timezone
import hashlib
from importlib import metadata
import io
import json
import os
from pathlib import Path
import platform
import re
import subprocess
from typing import Any
from urllib.parse import parse_qsl, urlsplit, urlunsplit
import zipfile


MAX_ARTIFACT_BYTES = 20 * 1024 * 1024
MAX_PACKAGE_BYTES = 200 * 1024 * 1024
_OMIT_DIRS = {'.git', '.venv', 'venv', 'env', 'node_modules', '__pycache__',
              '.pytest_cache', '.ipynb_checkpoints', 'reviews', 'review_packages'}
_SECRET_NAMES = re.compile(r'(^\.env(?:\.|$)|credentials?|cookies?|secrets?|tokens?|id_rsa|id_ed25519|\.pem$|\.key$)', re.I)
_SECRET_KEYS = re.compile(r'^(?:password|passwd|secret|token|access_token|refresh_token|api_key|apikey|authorization|cookie|credentials?)$', re.I)
_ACCESS_QUERY = re.compile(r'(token|secret|password|signature|credential|authorization|api.?key|^key$|^sig$|^code$|^auth$)', re.I)
_URL = re.compile(r'https?://[^\s<>\"\x27]+', re.I)
_TEXT_EXTENSIONS = {'.py', '.md', '.txt', '.json', '.csv', '.tsv', '.log', '.ipynb',
                    '.toml', '.yml', '.yaml', '.cxc', '.a3m', '.fasta', '.fa', '.faa',
                    '.pdb', '.cif', '.mmcif', '.svg'}
_IMAGE_EXTENSIONS = {'.png', '.jpg', '.jpeg'}
_SOURCE_DIRS = {'project1', 'tests', 'docs'}
_ROOT_FILES = {'README.md', 'START_HERE.md', 'AGENTS.md', 'requirements.txt',
               'requirements-dev.txt', 'pyproject.toml', 'pytest.ini',
               'Project1_Server_Workflow.ipynb', 'project1_config.json',
               'config.json', 'target_manifest.json', 'target_manifest.csv',
               'targets.json', 'manifest.json', '.gitignore'}
_QUESTIONS = (
    ('accuracy', "Did changing the MSA change the prediction's measured accuracy for each target?"),
    ('alignments', 'How did the two alignments differ where both are available?'),
    ('confidence_accuracy', 'Did confidence and measured accuracy move together?'),
    ('validation', 'Were controls, mapping, and metric implementations validated?'),
)
_ANSWER_STATES = {'ANSWERED_FOR_AVAILABLE_RUNS', 'PARTIALLY_ANSWERED', 'UNANSWERED'}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def redact_text(text: str) -> str:
    """Remove access-bearing URLs from package copies, retaining public sources.

    Query parameters are retained only if none resemble authorization material.
    URLs with embedded credentials or likely credential-bearing path segments
    are replaced in full. Originals on disk are never changed.
    """
    def clean(match: re.Match[str]) -> str:
        raw = match.group(0)
        try:
            parts = urlsplit(raw)
            path_access = re.search(r'/(?:token|access_token|auth|apikey|api_key|secret)/', parts.path, re.I)
            query_access = any(_ACCESS_QUERY.search(key) for key, _ in parse_qsl(parts.query))
            fragment_access = bool(_ACCESS_QUERY.search(parts.fragment)) if parts.fragment else False
            if parts.username or parts.password or path_access or query_access or fragment_access:
                return '[REDACTED_ACCESS_URL]'
            # A URL fragment is not useful for provenance and may hold an opaque token.
            return urlunsplit((parts.scheme, parts.netloc, parts.path, parts.query, ''))
        except ValueError:
            return '[REDACTED_INVALID_URL]'
    text = _URL.sub(clean, text)
    # Header/log forms. JSON objects are additionally sanitized structurally below.
    text = re.sub(r'(?im)^(\s*(?:authorization|cookie|set-cookie|password|api[_-]?key|access[_-]?token)\s*[:=]\s*).+$',
                  r'\1[REDACTED]', text)
    return text


def _sanitize(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): '[REDACTED]' if _SECRET_KEYS.fullmatch(str(key)) else _sanitize(item)
                for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_sanitize(item) for item in value]
    if isinstance(value, str):
        return redact_text(value)
    if isinstance(value, Path):
        return str(value)
    return value


def _json_bytes(value: Any) -> bytes:
    return (json.dumps(_sanitize(value), indent=2, ensure_ascii=False, default=str) + '\n').encode('utf-8')


def _is_within(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def _walk(root: Path):
    """Do not traverse environments, symlinks, credentials, or old packages."""
    if not root.is_dir():
        return
    for directory, names, files in os.walk(root, followlinks=False):
        names[:] = sorted(name for name in names if name.lower() not in _OMIT_DIRS
                          and not _SECRET_NAMES.search(name)
                          and not Path(directory, name).is_symlink())
        for name in sorted(files):
            path = Path(directory, name)
            if not path.is_symlink() and not _SECRET_NAMES.search(name):
                yield path


def _source_allowed(relative: Path) -> bool:
    if len(relative.parts) == 1:
        return relative.name in _ROOT_FILES or relative.suffix.lower() in {'.fasta', '.fa', '.faa'}
    if relative.parts[0] in _SOURCE_DIRS:
        # Include synthetic test SOURCE only, never generated fixture data.
        if relative.parts[0] == 'tests':
            return relative.suffix == '.py'
        return relative.suffix.lower() in _TEXT_EXTENSIONS | _IMAGE_EXTENSIONS
    if relative.parts[0] in {'config', 'configs'}:
        return relative.suffix.lower() in {'.json', '.csv', '.toml', '.yaml', '.yml'}
    if relative.parts[0] == 'outputs':
        if len(relative.parts) > 1 and relative.parts[1] == 'analyses':
            # Historical input-hash directories are retained locally but must
            # not look like current evidence in a newly generated package.
            return False
        return relative.suffix.lower() in _TEXT_EXTENSIONS | _IMAGE_EXTENSIONS
    return False


def _artifact_copy(path: Path) -> bytes:
    if path.suffix.lower() in {'.json', '.ipynb'}:
        try:
            value = json.loads(path.read_text(encoding='utf-8-sig'))
            if path.suffix.lower() == '.ipynb' and isinstance(value, dict):
                # Inputs reproduce the notebook; transient outputs can contain secrets
                # or stale data and are intentionally excluded from the review copy.
                for cell in value.get('cells', []):
                    if cell.get('cell_type') == 'code':
                        cell['outputs'] = []
                        cell['execution_count'] = None
            return _json_bytes(value)
        except (ValueError, UnicodeError):
            pass
    if path.suffix.lower() in _TEXT_EXTENSIONS or path.name == '.gitignore':
        return redact_text(path.read_text(encoding='utf-8-sig')).encode('utf-8')
    return path.read_bytes()


def _git_state(repo_root: Path) -> dict:
    state = {'revision': None, 'dirty': None, 'status': 'UNAVAILABLE'}
    try:
        revision = subprocess.run(['git', '-C', str(repo_root), 'rev-parse', 'HEAD'],
                                  capture_output=True, text=True, timeout=10, check=True)
        dirty = subprocess.run(['git', '-C', str(repo_root), 'status', '--porcelain', '--untracked-files=normal'],
                               capture_output=True, text=True, timeout=10, check=True)
        state.update(revision=revision.stdout.strip(), dirty=bool(dirty.stdout.strip()),
                     status='RECORDED', changes=dirty.stdout.splitlines())
    except (OSError, subprocess.SubprocessError):
        state['reason'] = 'Git revision or working-tree status could not be read.'
    return state


def _versions() -> dict:
    versions = {'python': platform.python_version(), 'platform': platform.platform()}
    for package in ('numpy', 'pandas', 'biopython', 'matplotlib', 'pytest', 'nbformat', 'nbclient', 'ipykernel', 'jupyterlab', 'tmtools'):
        try:
            versions[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            versions[package] = 'NOT_INSTALLED'
    return versions


def _target_rows(status: dict) -> list[dict]:
    targets = status.get('targets', {})
    if isinstance(targets, dict):
        return [dict(row, target=key) if isinstance(row, dict) else {'target': key, 'state': row}
                for key, row in targets.items()]
    if isinstance(targets, list):
        return [row for row in targets if isinstance(row, dict)]
    return []


def _real_summary(status: dict) -> dict:
    """Omit rows explicitly labeled synthetic; never collect test directories."""
    copied = dict(status)
    copied['targets'] = [row for row in _target_rows(status)
                         if not row.get('synthetic', False)
                         and str(row.get('data_kind', '')).lower() not in {'synthetic', 'fixture', 'test'}]
    return copied


def _evidence_exists(evidence: str, repo_root: Path, data_root: Path) -> bool:
    path = Path(evidence)
    candidates = [path] if path.is_absolute() else [repo_root / path, data_root / path]
    return any(candidate.is_file() and not candidate.is_symlink()
               and (_is_within(candidate.resolve(), repo_root) or _is_within(candidate.resolve(), data_root))
               and not any(part.lower() in {'tests', 'fixtures', 'synthetic'} for part in candidate.parts)
               for candidate in candidates)


def _review_text(status: dict, repo_root: Path, data_root: Path, created: str) -> str:
    lines = ['# Project 1 evidence review', '', f'Created (UTC): {created}', '',
             'This package records the available local evidence. Prepared files are not proof of server submission.',
             'Synthetic tests validate software behavior and are not measured study results.', '',
             'Status snapshot: `provenance/status_snapshot.json`. Artifact hashes and exclusions: '
             '`provenance/artifact_inventory.csv`. Software and Git state: `provenance/environment.json`.', '']
    supplied = status.get('review_answers', {})
    for number, (key, question) in enumerate(_QUESTIONS, 1):
        candidate = supplied.get(key, {}) if isinstance(supplied, dict) else {}
        evidence = candidate.get('evidence', []) if isinstance(candidate, dict) else []
        if isinstance(evidence, str):
            evidence = [evidence]
        evidence = [str(item) for item in evidence]
        valid = [item for item in evidence if _evidence_exists(item, repo_root, data_root)]
        state = candidate.get('status', 'UNANSWERED') if isinstance(candidate, dict) else 'UNANSWERED'
        if state not in _ANSWER_STATES or (state != 'UNANSWERED' and not valid):
            state = 'UNANSWERED'
        lines.extend([f'## {number}. {question}', '', f'**{state}**', ''])
        if state != 'UNANSWERED':
            lines.append(str(candidate.get('answer', 'Evidence is available; see the cited files for target-specific results.')))
        elif key == 'accuracy':
            lines.append('No verified paired accuracy result is established by this review. Both primary models, '
                         'verified matched controls, a confirmed experimental reference, and a common residue mask are required.')
        elif key == 'alignments':
            lines.append('No verified comparison of both real alignments is established. An unavailable returned AUTO MSA '
                         'must remain unavailable; the target sequence alone is not a substitute alignment.')
        elif key == 'confidence_accuracy':
            lines.append('No relationship is established without paired real confidence outputs and measured reference accuracy. '
                         'Confidence measures do not establish accuracy on their own.')
        else:
            lines.append('Real-run controls and residue mappings remain unverified unless explicitly supported by the imported evidence. '
                         'See `repo/outputs/tests.log` when present for synthetic software checks; such tests do not validate a server run.')
        if valid:
            lines.append('Evidence: ' + ', '.join(f'`{item}`' for item in valid) + '.')
        invalid = [item for item in evidence if item not in valid]
        if invalid:
            lines.append('Unverified or unavailable evidence paths: ' + ', '.join(f'`{item}`' for item in invalid) + '.')
        if isinstance(candidate, dict) and candidate.get('limitations'):
            lines.append('Limitations: ' + str(candidate['limitations']))
        lines.append('')
    lines.extend(['## Metric availability', '',
                  'Use the generated availability tables in `repo/outputs/` and '
                  '`provenance/metric_availability.csv`. Missing values are missing, never zero. '
                  'A planned metric lacking both a real value and evidence remains unavailable.', '',
                  '## Interpretation limits', '',
                  'Returned models are not independent protein replicates. No p-values across such models are justified here. '
                  'Three targets do not establish general superiority, novelty, or a recreated blind CASP benchmark. '
                  'Disabling structural templates does not prove training-data independence.', '',
                  '## Package policy', '',
                  'The package uses an allowlist, strips notebook outputs, and redacts recognized access-bearing URLs '
                  'and credential fields in text copies. Original files are unchanged. Source tests may contain synthetic examples; '
                  'test fixtures are not study observations. Raw downloads are inventoried rather than automatically copied. '
                  'Oversized artifacts are omitted intact and listed with full local location, size, and SHA-256. '
                  'No alignment is truncated. Virtual environments, Git internals, credential files, and prior review ZIPs are excluded.', ''])
    return redact_text('\n'.join(lines))


def _metric_rows(status: dict, repo_root: Path, data_root: Path) -> list[dict]:
    """Complete planned coverage even when no analysis has yet been possible."""
    metrics = {
        'alignment_coverage': 'Validated real alignment with query-coordinate mapping',
        'identity_to_query': 'Validated real alignment',
        'sequence_redundancy': 'Validated real alignment',
        'normalized_BLOSUM62_sum_of_pairs': 'Validated real alignment; documented eligible-pair denominator',
        'Neff': 'Declared identity/coverage definition and bounded exact or labeled approximate implementation',
        'residue_mean_pLDDT': 'Real confidence output with atom-to-residue mapping',
        'pTM': 'Real server confidence output',
        'PAE': 'Real server confidence output and residue mapping',
        'CA_RMSD': 'Verified paired primary models, confirmed reference, frozen common mask',
        'TM_score': 'Standard optimized TM-score tool and verified common reference normalization',
        'lDDT_CA': 'Documented lDDT-C-alpha implementation and frozen common mask',
        'GDT_TS': 'Standard GDT-TS tool with repeated fitting',
        'MolProbity': 'Supported installed MolProbity tool; separate geometry validation',
    }
    declared = status.get('metric_availability', [])
    declared = declared if isinstance(declared, list) else []
    lookup = {(str(row.get('target')), str(row.get('metric'))): row
              for row in declared if isinstance(row, dict) and not row.get('synthetic', False)}
    result = []
    targets = _target_rows(_real_summary(status))
    for target in targets:
        target_id = target.get('target', target.get('target_id', target.get('id', 'UNRESOLVED')))
        for metric, requirement in metrics.items():
            supplied = lookup.get((str(target_id), metric), {})
            state = supplied.get('status', 'UNAVAILABLE')
            evidence = str(supplied.get('evidence', ''))
            if state in {'AVAILABLE', 'VALIDATED', 'COMPUTED'} and not _evidence_exists(evidence, repo_root, data_root):
                state = 'UNAVAILABLE'
            result.append({'target': target_id, 'metric': metric, 'status': state,
                           'requirement': supplied.get('requirement', requirement),
                           'evidence': evidence,
                           'note': supplied.get('note', 'No verified real-data result supplied to the review generator.')})
    return result


def _csv_bytes(rows: list[dict], fields: list[str]) -> bytes:
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)
    return redact_text(stream.getvalue()).encode('utf-8')


def create_review(repo_root: str | Path, data_root: str | Path, status: dict,
                  output_root: str | Path | None = None) -> Path:
    """Write a timestamped review ZIP without modifying study originals.

    ``review_artifacts`` is an optional list of files/directories under repo/data roots that
    are approved for inclusion (validated alignments, primary structures,
    references, mappings, confidence inputs). Files above 20 MiB and packages
    above 200 MiB are inventoried instead, without truncation. The optional
    ``metric_availability`` may supply rows with target, metric, status, evidence,
    requirement, and note. Claimed available metrics require an existing real
    evidence file. ``review_answers`` uses keys accuracy, alignments,
    confidence_accuracy, validation, each with status, answer, evidence, and
    limitations. Non-UNANSWERED claims require existing nonsynthetic evidence.
    """
    repo_root, data_root = Path(repo_root).resolve(), Path(data_root).resolve()
    output_root = Path(output_root).resolve() if output_root else repo_root / 'outputs' / 'reviews'
    output_root.mkdir(parents=True, exist_ok=True)
    now = datetime.now(timezone.utc)
    archive_path = output_root / f'project1_review_{now.strftime("%Y%m%dT%H%M%S%fZ")}.zip'
    real_status = _real_summary(status)
    # Keep the current review readable without unpacking a ZIP. Serialize once
    # so the on-disk snapshots and ZIP provenance describe the same run.
    review_content = _review_text(real_status, repo_root, data_root, now.isoformat()).encode('utf-8')
    versions, git_state = _versions(), _git_state(repo_root)
    metric_content = _csv_bytes(_metric_rows(real_status, repo_root, data_root),
                                ['target', 'metric', 'status', 'requirement', 'evidence', 'note'])
    generated = repo_root / 'outputs'
    generated.mkdir(parents=True, exist_ok=True)
    for name, content in {
        'REVIEW.md': review_content,
        'dependency_versions.json': _json_bytes(versions),
        'git_provenance.json': _json_bytes(git_state),
        'metric_availability.csv': metric_content,
    }.items():
        (generated / name).write_bytes(content)
    explicit, explicit_dirs = set(), set()
    artifact_paths = list(status.get('review_artifacts', []))
    if 'review_artifacts' not in status:
        current_path = repo_root / 'outputs' / 'current.json'
        if current_path.is_file():
            try:
                artifact_paths.extend(json.loads(current_path.read_text(encoding='utf-8')).get('artifacts', []))
            except (ValueError, UnicodeError, OSError):
                pass
    # Include only currently prepared submission inputs, not stale requests
    # left behind when a sequence or alignment was replaced.
    for target in _target_rows(real_status):
        for stage in ('AUTO', 'DEEPMSA', 'deepmsa'):
            row = target.get(stage, {})
            if isinstance(row, dict):
                artifact_paths.extend(row[key] for key in ('request_path', 'fasta_path', 'sequence_path')
                                      if row.get(key))
    for item in artifact_paths:
        path = Path(item)
        path = (repo_root / path).resolve() if not path.is_absolute() else path.resolve()
        if _is_within(path, repo_root) or _is_within(path, data_root):
            if path.is_dir():
                explicit_dirs.add(path)
            else:
                explicit.add(path)
    inventory, total, seen = [], 0, set()
    roots = [(repo_root, 'repo')]
    if data_root != repo_root and not _is_within(data_root, repo_root):
        roots.append((data_root, 'data'))
    with zipfile.ZipFile(archive_path, mode='x', compression=zipfile.ZIP_DEFLATED) as package:
        for root, prefix in roots:
            for path in _walk(root):
                resolved = path.resolve()
                if resolved in seen or _is_within(resolved, output_root):
                    continue
                seen.add(resolved)
                relative = path.relative_to(root)
                if any(part.lower() in {'fixtures', 'synthetic'} for part in relative.parts):
                    continue
                source = prefix == 'repo' and _source_allowed(relative)
                safe_extension = path.suffix.lower() in _TEXT_EXTENSIONS | _IMAGE_EXTENSIONS
                explicitly_selected = resolved in explicit or any(_is_within(resolved, directory) for directory in explicit_dirs)
                allowed = source or (explicitly_selected and safe_extension)
                size = path.stat().st_size
                if path.name.startswith('project1_review_') and path.suffix == '.zip':
                    continue
                reason = '' if allowed else 'outside package allowlist (original retained)'
                if size > MAX_ARTIFACT_BYTES:
                    reason = 'exceeds per-artifact size limit (original retained intact)'
                elif allowed and total + size > MAX_PACKAGE_BYTES:
                    reason = 'exceeds package size budget (original retained intact)'
                # Every allowed item and every retained original is accounted for;
                # secret/system files were excluded before opening them.
                record = {'location': str(resolved), 'size_bytes': size,
                          'source_sha256': _sha256(path), 'archive_path': '',
                          'packaged_sha256': '', 'state': 'EXCLUDED', 'reason': reason}
                if not reason:
                    try:
                        content = _artifact_copy(path)
                    except UnicodeError:
                        record['reason'] = 'text decoding failed; original retained'
                    else:
                        if len(content) > MAX_ARTIFACT_BYTES or total + len(content) > MAX_PACKAGE_BYTES:
                            record['reason'] = 'sanitized copy exceeds package/artifact size limit (original retained intact)'
                        else:
                            destination = f'{prefix}/{relative.as_posix()}'
                            package.writestr(destination, content)
                            record.update(archive_path=destination, packaged_sha256=hashlib.sha256(content).hexdigest(),
                                          state='INCLUDED', reason='sanitized copy; compare source and packaged hashes')
                            total += len(content)
                inventory.append(record)
        package.writestr('REVIEW.md', review_content)
        package.writestr('provenance/status_snapshot.json', _json_bytes(real_status))
        package.writestr('provenance/environment.json', _json_bytes({'created_utc': now.isoformat(),
                          'dependency_versions': versions, 'git': git_state}))
        package.writestr('provenance/artifact_inventory.csv', _csv_bytes(inventory,
                          ['location', 'size_bytes', 'source_sha256', 'archive_path', 'packaged_sha256', 'state', 'reason']))
        package.writestr('provenance/metric_availability.csv', metric_content)
    return archive_path
