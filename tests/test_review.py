"""Synthetic review-package tests. No files here are study results."""
import csv
import io
import json
import zipfile

from project1 import review


def _write(root, relative, value):
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value, encoding='utf-8')
    return path


def test_review_safelists_and_redacts_copies_without_touching_original(tmp_path):
    repo, data = tmp_path / 'repo', tmp_path / 'data'
    source = _write(repo, 'outputs/provenance.json', json.dumps({
        'url': 'https://example.org/download?token=topsecret', 'password': 'also_secret',
        'public': 'https://example.org/documentation'}))
    original = source.read_bytes()
    _write(repo, '.env', 'secret=hide')
    _write(repo, '.git/config', 'password=hide')
    _write(repo, '.venv/secret.txt', 'hide')
    _write(repo, 'tests/test_metrics.py', '# synthetic software tests\n')
    _write(repo, 'tests/fixtures/fake_result.csv', 'fake,42\n')
    raw = _write(data, 'downloads/original.a3m', '>synthetic_query\nACD\n')
    archive = review.create_review(repo, data, {'targets': {'REAL_TARGET': {'state': 'WAITING'}}})
    assert source.read_bytes() == original
    with zipfile.ZipFile(archive) as package:
        names = package.namelist()
        assert 'repo/tests/test_metrics.py' in names
        assert not any('.env' in name or '.git/' in name or '.venv/' in name or 'fixtures/' in name for name in names)
        copied = json.loads(package.read('repo/outputs/provenance.json'))
        assert copied['url'] == '[REDACTED_ACCESS_URL]'
        assert copied['password'] == '[REDACTED]'
        assert copied['public'] == 'https://example.org/documentation'
        report = package.read('REVIEW.md').decode()
        assert report.count('**UNANSWERED**') == 4
        assert 'topsecret' not in report
        rows = list(csv.DictReader(io.StringIO(package.read('provenance/artifact_inventory.csv').decode())))
        entry = next(row for row in rows if row['location'] == str(raw.resolve()))
        assert entry['state'] == 'EXCLUDED'
        assert entry['source_sha256'] == review._sha256(raw)


def test_large_alignment_is_inventoried_not_truncated(tmp_path, monkeypatch):
    repo, data = tmp_path / 'repo', tmp_path / 'data'
    repo.mkdir()
    alignment = _write(data, 'validated/real.a3m', '>synthetic\n' + 'A' * 100)
    monkeypatch.setattr(review, 'MAX_ARTIFACT_BYTES', 20)
    archive = review.create_review(repo, data, {'review_artifacts': [str(alignment)]})
    with zipfile.ZipFile(archive) as package:
        assert not any(name.endswith('.a3m') for name in package.namelist())
        rows = list(csv.DictReader(io.StringIO(package.read('provenance/artifact_inventory.csv').decode())))
        entry = next(row for row in rows if row['location'] == str(alignment.resolve()))
        assert entry['size_bytes'] == str(alignment.stat().st_size)
        assert 'per-artifact size' in entry['reason']


def test_reviews_unique_and_do_not_repackage_prior_reviews(tmp_path):
    repo = tmp_path / 'repo'
    _write(repo, 'README.md', 'Documentation\n')
    first = review.create_review(repo, repo / 'data', {})
    second = review.create_review(repo, repo / 'data', {})
    assert first != second
    with zipfile.ZipFile(second) as package:
        assert not any(name.endswith('.zip') for name in package.namelist())


def test_synthetic_status_excluded_and_unproven_answer_downgraded(tmp_path):
    repo = tmp_path / 'repo'
    _write(repo, 'README.md', 'Test\n')
    status = {'targets': {'REAL': {'state': 'WAITING'}, 'FAKE': {'synthetic': True, 'rmsd': 12345}},
              'review_answers': {'accuracy': {'status': 'ANSWERED_FOR_AVAILABLE_RUNS',
                                             'answer': 'UNSUPPORTED CONCLUSION', 'evidence': ['missing.csv']}}}
    archive = review.create_review(repo, repo / 'data', status)
    with zipfile.ZipFile(archive) as package:
        snapshot = package.read('provenance/status_snapshot.json').decode()
        assert 'FAKE' not in snapshot and '12345' not in snapshot
        report = package.read('REVIEW.md').decode()
        assert 'UNSUPPORTED CONCLUSION' not in report
        assert '**UNANSWERED**' in report


def test_notebook_copy_drops_transient_outputs(tmp_path):
    repo = tmp_path / 'repo'
    notebook = {'cells': [{'cell_type': 'code', 'source': ['print(1)'], 'execution_count': 4,
                           'outputs': [{'text': ['private transient output']}]}], 'nbformat': 4}
    source = _write(repo, 'Project1_Server_Workflow.ipynb', json.dumps(notebook))
    archive = review.create_review(repo, repo / 'data', {})
    with zipfile.ZipFile(archive) as package:
        copied = json.loads(package.read('repo/Project1_Server_Workflow.ipynb'))
        assert copied['cells'][0]['outputs'] == []
        assert copied['cells'][0]['execution_count'] is None
    assert json.loads(source.read_text()) == notebook


def test_access_urls_and_authorization_redacted():
    inputs = ['https://user:pass@example.org/x', 'https://example.org/x?X-Amz-Signature=secret',
              'https://example.org/token/secret', 'https://example.org/x?sig=secret']
    for url in inputs:
        assert review.redact_text(url) == '[REDACTED_ACCESS_URL]'
    assert 'abc123' not in review.redact_text('Authorization: Bearer abc123\n')


def test_metric_availability_requires_existing_nonsynthetic_evidence(tmp_path):
    repo = tmp_path / 'repo'
    _write(repo, 'outputs/measured.csv', 'target,rmsd\nREAL,1.2\n')
    _write(repo, 'tests/fixtures/fake.csv', 'fake,42\n')
    status = {'targets': {'REAL': {}}, 'metric_availability': [
        {'target': 'REAL', 'metric': 'CA_RMSD', 'status': 'AVAILABLE', 'evidence': 'outputs/measured.csv'},
        {'target': 'REAL', 'metric': 'TM_score', 'status': 'AVAILABLE', 'evidence': 'tests/fixtures/fake.csv'},
        {'target': 'REAL', 'metric': 'pTM', 'status': 'AVAILABLE', 'evidence': 'missing.csv'}]}
    archive = review.create_review(repo, repo / 'data', status)
    with zipfile.ZipFile(archive) as package:
        rows = {row['metric']: row for row in csv.DictReader(io.StringIO(package.read('provenance/metric_availability.csv').decode()))}
        assert rows['CA_RMSD']['status'] == 'AVAILABLE'
        assert rows['TM_score']['status'] == 'UNAVAILABLE'
        assert rows['pTM']['status'] == 'UNAVAILABLE'
        assert rows['Neff']['status'] == 'UNAVAILABLE'


def test_current_artifact_directories_include_nested_files_and_exclude_stale_analyses(tmp_path):
    repo, data = tmp_path / 'repo', tmp_path / 'data'
    current = _write(repo, 'outputs/analyses/REAL/current/alignment/metrics.csv', 'data\n')
    _write(repo, 'outputs/analyses/REAL/old/alignment/metrics.csv', 'stale\n')
    external = _write(data, 'validated/REAL/current/model.cif', '# synthetic test\n')
    stale_request = _write(repo, 'submissions/REAL/old_request.json', '{}')
    request = _write(repo, 'submissions/REAL/current_request.json', '{}')
    status = {'targets': {'REAL': {'AUTO': {'request_path': str(request)}}},
              'review_artifacts': [str(current.parent.parent), str(external.parent)]}
    archive = review.create_review(repo, data, status)
    with zipfile.ZipFile(archive) as package:
        names = package.namelist()
        assert 'repo/outputs/analyses/REAL/current/alignment/metrics.csv' in names
        assert 'data/validated/REAL/current/model.cif' in names
        assert 'repo/submissions/REAL/current_request.json' in names
        assert not any('/old/' in name for name in names)
        assert not any(name.endswith(stale_request.name) for name in names)


def test_explicit_empty_artifact_list_overrides_old_current_snapshot(tmp_path):
    repo = tmp_path / 'repo'
    stale = _write(repo, 'outputs/analyses/REAL/old/metrics.csv', 'stale\n')
    _write(repo, 'outputs/current.json', json.dumps({'artifacts': [str(stale.parent)]}))
    archive = review.create_review(repo, repo / 'data', {'review_artifacts': []})
    with zipfile.ZipFile(archive) as package:
        assert not any('/analyses/' in name for name in package.namelist())


def test_current_review_and_provenance_written_before_packaging(tmp_path):
    repo = tmp_path / 'repo'
    archive = review.create_review(repo, repo / 'data', {'targets': {'REAL': {}}})
    with zipfile.ZipFile(archive) as package:
        assert (repo / 'outputs/REVIEW.md').read_bytes() == package.read('REVIEW.md')
        assert package.read('repo/outputs/REVIEW.md') == package.read('REVIEW.md')
        assert (repo / 'outputs/metric_availability.csv').read_bytes() == package.read('provenance/metric_availability.csv')
        environment = json.loads(package.read('provenance/environment.json'))
        versions = json.loads((repo / 'outputs/dependency_versions.json').read_text())
        git_state = json.loads((repo / 'outputs/git_provenance.json').read_text())
        assert versions == environment['dependency_versions']
        assert git_state == environment['git']
        assert json.loads(package.read('repo/outputs/dependency_versions.json')) == versions
        assert json.loads(package.read('repo/outputs/git_provenance.json')) == git_state
