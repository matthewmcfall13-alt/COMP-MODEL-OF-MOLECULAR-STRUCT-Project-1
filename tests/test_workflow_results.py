"""End-to-end local orchestration with SYNTHETIC files isolated under tmp_path.

No server access, real target identifiers, downloads, or study directories used.
"""
import json
from pathlib import Path
import zipfile

import numpy as np
from Bio.PDB import MMCIFIO, PDBParser
import pytest

from project1.inputs import file_hash
from project1.metrics import _write_ca_pdb
from project1.workflow import run


TARGET = "SYNTHETIC_TEST_TARGET"
SEQUENCE = "ACDEF"
A3M = ">SYNTHETIC_QUERY\nACDEF\n>UniRef90_SYNTHETIC_HIT\nACdD-F\n"


@pytest.fixture
def study(tmp_path):
    root = tmp_path / "synthetic_repo"
    root.mkdir()
    fasta = root / "source.fasta"
    fasta.write_bytes(b">SYNTHETIC_TEST_TARGET original synthetic header\r\nACDEF\r\n")
    config = {"repo_root": ".", "data_root": str(tmp_path / "external synthetic data with spaces"),
              "selected_targets": [TARGET], "pilot": TARGET, "seed": 42,
              "input_overrides": {}, "tmscore_executable": "__absent_synthetic_TMscore__",
              "targets": {TARGET: {"deepmsa": {}, "AUTO": {}, "DEEPMSA": {}, "reference": {}}}}
    path = root / "config.json"
    path.write_text(json.dumps(config))
    return root, config, path


def save_config(study):
    study[2].write_text(json.dumps(study[1]))


def initial_with_a3m(study):
    state = run(study[2], package=False)
    folder = Path(state["targets"][TARGET]["deepmsa"]["download_dir"])
    alignment = folder / "final_SYNTHETIC.a3m"
    alignment.write_text(A3M)
    return alignment, run(study[2], package=False)


def synthetic_public_zip(folder, prepared_request, tmp_path, wrong_structure=False):
    """Produce a public-server-shaped archive, explicitly synthetic content."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    coordinates = np.array([[0, 0, 0], [3, 0, 0], [3, 4, 0], [3, 4, 5], [1, 2, 6]], dtype=float)
    pdb = tmp_path / "model_SYNTHETIC.pdb"
    _write_ca_pdb(pdb, list(range(1, 6)), "CCDEF" if wrong_structure else SEQUENCE, coordinates)
    cif = tmp_path / "model_SYNTHETIC.cif"
    writer = MMCIFIO()
    writer.set_structure(PDBParser(QUIET=True).get_structure("SYNTHETIC", pdb))
    writer.save(str(cif))
    archive_path = folder / "SYNTHETIC_PUBLIC_SERVER_RESULT.zip"
    with zipfile.ZipFile(archive_path, "w") as archive:
        archive.writestr("nested/fold_SYNTHETIC_job_request.json", json.dumps(prepared_request))
        archive.writestr("nested/fold_SYNTHETIC_model_0.cif", cif.read_text())
        archive.writestr("nested/fold_SYNTHETIC_summary_confidences_0.json", json.dumps({"ranking_score": 0.9, "ptm": 0.8}))
        archive.writestr("nested/fold_SYNTHETIC_full_data_0.json", json.dumps({
            "atom_plddts": [80, 81, 82, 83, 84], "atom_chain_ids": ["A"] * 5,
            "pae": [[0 if i == j else 2 for j in range(5)] for i in range(5)],
            "token_chain_ids": ["A"] * 5, "token_res_ids": [1, 2, 3, 4, 5]}))
    return archive_path, pdb


def create_valid_pair(study, tmp_path):
    _, state = initial_with_a3m(study)
    target = state["targets"][TARGET]
    for condition in ("AUTO", "DEEPMSA"):
        expected = json.loads(Path(target[condition]["request_path"]).read_text())
        _, reference = synthetic_public_zip(target[condition]["download_dir"], expected, tmp_path)
    study[1]["targets"][TARGET]["reference"] = {
        "path": str(reference), "chain": "A", "confirmed": True, "experimental": True,
        "assignment_evidence": "SYNTHETIC TEST ONLY: explicit temporary fixture assignment, no scientific claim"}
    save_config(study)
    return run(study[2], package=False)


def test_a3m_progression_preserves_full_insertions_and_request_idempotency(study):
    original_hash = file_hash(study[0] / "source.fasta")
    alignment, first = initial_with_a3m(study)
    row = first["targets"][TARGET]
    assert row["deepmsa"]["state"] == "VALIDATED"
    assert row["DEEPMSA"]["state"] == "PREPARED"
    deep_path = Path(row["DEEPMSA"]["request_path"])
    content, timestamp = deep_path.read_bytes(), deep_path.stat().st_mtime_ns
    request = json.loads(content)
    assert request[0]["sequences"][0]["proteinChain"]["unpairedMsa"] == alignment.read_bytes().decode("utf-8")
    analysis_copy = Path(row["deepmsa"]["analysis_directory"]) / "query_coordinates.fasta"
    assert "ACD-F" in analysis_copy.read_text()
    assert alignment.read_text() == A3M
    second = run(study[2], package=False)
    assert second["targets"][TARGET]["DEEPMSA"]["request_path"] == str(deep_path)
    assert deep_path.read_bytes() == content
    assert deep_path.stat().st_mtime_ns == timestamp
    assert file_hash(study[0] / "source.fasta") == original_hash
    assert len(list((study[0] / "submissions" / TARGET).glob("*.json"))) == 2
    assert list(json.loads((study[0] / "target_manifest.json").read_text())["targets"]) == [TARGET]


def test_stage_directories_are_distinct_on_case_insensitive_filesystems(study):
    row = run(study[2], package=False)["targets"][TARGET]
    directories = [row[stage]["download_dir"].casefold() for stage in ("deepmsa", "AUTO", "DEEPMSA")]
    assert len(set(directories)) == 3, "DeepMSA2 and AF custom-MSA downloads cannot differ by case only on Windows"


def test_external_input_override_and_manual_search_metadata(study, tmp_path):
    external_fasta = tmp_path / "source with spaces.txt"
    external_fasta.write_text(f">{TARGET} custom synthetic metadata\n{SEQUENCE}\n")
    external_msa = tmp_path / "unrelated_filename.a3m"
    external_msa.write_text(A3M)
    study[1]["input_overrides"] = {TARGET: str(external_fasta)}
    study[1]["targets"][TARGET]["deepmsa"] = {
        "input_path": str(external_msa), "actual_mode": "medium",
        "job_id": "SYNTHETIC_ID", "submitted_at": "2026-09-22T00:00:00Z",
        "search_metadata": {"database_release": "SYNTHETIC_RECORDED_METADATA"}}
    save_config(study)
    result = run(study[2], package=False)["targets"][TARGET]
    manifest = json.loads((study[0] / "target_manifest.json").read_text())["targets"][TARGET]
    assert manifest["source_path"] == str(external_fasta.resolve())
    assert result["deepmsa"]["state"] == "VALIDATED"
    assert result["deepmsa"]["actual_mode"] == "medium"
    assert result["deepmsa"]["manual"]["job_id"] == "SYNTHETIC_ID"
    assert result["deepmsa"]["search_metadata"]["database_release"] == "SYNTHETIC_RECORDED_METADATA"


def test_custom_msa_submission_metadata_survives_rerun(study):
    _, state = initial_with_a3m(study)
    target_options = study[1]["targets"][TARGET]
    target_options["DEEPMSA"] = {
        "submitted_at": "2026-09-22T01:02:03Z", "job_id": "SYNTHETIC_MANUAL_JOB",
        "submitted_request_sha256": state["targets"][TARGET]["DEEPMSA"]["request_sha256"]}
    save_config(study)
    first = run(study[2], package=False)["targets"][TARGET]["DEEPMSA"]
    assert first["state"] == "SUBMITTED"
    target_options["DEEPMSA"] = {}
    save_config(study)
    second = run(study[2], package=False)["targets"][TARGET]["DEEPMSA"]
    assert second["state"] == "SUBMITTED"
    assert second["manual"]["job_id"] == "SYNTHETIC_MANUAL_JOB"


def test_explicit_submitted_hash_is_preserved_when_alignment_changes(study):
    alignment, state = initial_with_a3m(study)
    submitted_hash = state["targets"][TARGET]["DEEPMSA"]["request_sha256"]
    study[1]["targets"][TARGET]["DEEPMSA"] = {
        "submitted_at": "2026-09-22T01:02:03Z", "submitted_request_sha256": submitted_hash}
    save_config(study)
    assert run(study[2], package=False)["targets"][TARGET]["DEEPMSA"]["state"] == "SUBMITTED"
    study[1]["targets"][TARGET]["DEEPMSA"] = {}
    save_config(study)
    alignment.write_text(A3M.replace("ACdD-F", "ACeDEF"))
    changed = run(study[2], package=False)["targets"][TARGET]["DEEPMSA"]
    assert changed["state"] == "PREPARED"
    assert changed["manual"]["submitted_request_sha256"] == submitted_hash
    assert "differs" in changed["submission_warning"]


def test_changed_a3m_replaces_current_analysis_and_removed_download_waits(study):
    alignment, state = initial_with_a3m(study)
    old = state["targets"][TARGET]["deepmsa"]["analysis_directory"]
    alignment.write_text(A3M.replace("ACdD-F", "ACeDEF"))
    changed = run(study[2], package=False)
    new = changed["targets"][TARGET]["deepmsa"]["analysis_directory"]
    assert new != old
    assert Path(old).exists()  # Historical analysis preserved, never reused as current.
    assert old not in json.loads((study[0] / "outputs" / "current.json").read_text())["artifacts"]
    alignment.unlink()
    missing = run(study[2], package=False)
    assert missing["targets"][TARGET]["DEEPMSA"]["state"] == "WAITING_FOR_DEEPMSA"
    assert missing["targets"][TARGET]["deepmsa"]["state"] == "WAITING_FOR_DEEPMSA"
    current = json.loads((study[0] / "outputs" / "current.json").read_text())["artifacts"]
    assert new not in current and old not in current
    assert "request_path" not in missing["targets"][TARGET]["DEEPMSA"]


def test_multiple_a3m_candidates_block_until_explicit_selection(study):
    first, initial = initial_with_a3m(study)
    (first.parent / "another_SYNTHETIC.a3m").write_text(A3M.replace("ACdD-F", "ACeDEF"))
    blocked = run(study[2], package=False)["targets"][TARGET]
    assert blocked["deepmsa"]["state"] == "BLOCKED"
    assert blocked["DEEPMSA"]["state"] == "BLOCKED"
    assert len(blocked["deepmsa"]["candidate_inventory"]) == 2
    study[1]["targets"][TARGET]["deepmsa"]["a3m_selection"] = first.name
    save_config(study)
    selected = run(study[2], package=False)["targets"][TARGET]
    assert selected["deepmsa"]["state"] == "VALIDATED"
    assert "Explicit" in selected["deepmsa"]["selection_rationale"]


def test_valid_public_pair_reaches_accuracy_without_auto_alignment(study, tmp_path):
    state = create_valid_pair(study, tmp_path)
    row = state["targets"][TARGET]
    assert row["AUTO"]["state"] == "VALIDATED"
    assert row["DEEPMSA"]["state"] == "VALIDATED"
    assert row["accuracy"]["status"] == "ANALYZED"
    assert row["accuracy"]["scored_residue_count"] == 5
    assert row["AUTO"]["alignment_analysis"]["status"] == "UNAVAILABLE"
    assert row["accuracy"]["metrics"]["AUTO"]["rmsd_angstrom"] < 1e-8
    assert row["accuracy"]["metrics"]["AUTO"]["tm_score"] is None
    assert any(Path(item).name == "accuracy" for item in state["review_artifacts"])


def test_wrong_target_import_blocks_and_removes_prior_accuracy_from_current(study, tmp_path):
    state = create_valid_pair(study, tmp_path)
    assert state["targets"][TARGET]["accuracy"]["status"] == "ANALYZED"
    previous_accuracy = [item for item in state["review_artifacts"] if Path(item).name == "accuracy"]
    row = state["targets"][TARGET]
    expected = json.loads(Path(row["AUTO"]["request_path"]).read_text())
    synthetic_public_zip(row["AUTO"]["download_dir"], expected, tmp_path, wrong_structure=True)
    blocked = run(study[2], package=False)
    assert blocked["targets"][TARGET]["AUTO"]["state"] == "BLOCKED"
    assert blocked["targets"][TARGET]["accuracy"]["status"] != "ANALYZED"
    assert not blocked["targets"][TARGET]["AUTO"]["import"]["comparison_eligible"]
    current = json.loads((study[0] / "outputs" / "current.json").read_text())["artifacts"]
    assert all(path not in current for path in previous_accuracy)
    assert all(Path(path).exists() for path in previous_accuracy)
