"""SYNTHETIC ONLY: downloads generated in pytest tmp_path, never real results."""
import gzip
import io
import json
import stat
import tarfile
import zipfile
from pathlib import Path

import pytest

from project1.importers import ingest_files, import_alphafold, parse_confidence, verify_job_request


def request(sequence="AC", **changes):
    job = {"name": "SYNTHETIC_TEST_ONLY", "modelSeeds": ["42"], "sequences": [
        {"proteinChain": {"sequence": sequence, "count": 1, "useStructureTemplate": False}}
    ], "dialect": "alphafoldserver", "version": 1}
    job.update(changes)
    return [job]


def cif_text(names=("ALA", "CYS"), chain="A", values=(20, 20, 80), second_chain=None):
    # Two atoms in residue 1, one in residue 2 distinguish atom/residue means.
    rows = []
    for index, (resid, atom_name) in enumerate(((1, "N"), (1, "CA"), (2, "CA"))):
        use_chain = second_chain if resid == 2 and second_chain else chain
        rows.append(f"ATOM {index+1} C {atom_name} . {names[resid-1]} {use_chain} 1 {resid} ? {index}.0 0.0 0.0 1.0 {values[index]} {resid+9} {names[resid-1]} {use_chain} {atom_name} 1")
    return "data_SYNTHETIC_ONLY\nloop_\n" + "\n".join([
        "_atom_site.group_PDB", "_atom_site.id", "_atom_site.type_symbol", "_atom_site.label_atom_id",
        "_atom_site.label_alt_id", "_atom_site.label_comp_id", "_atom_site.label_asym_id",
        "_atom_site.label_entity_id", "_atom_site.label_seq_id", "_atom_site.pdbx_PDB_ins_code",
        "_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z", "_atom_site.occupancy",
        "_atom_site.B_iso_or_equiv", "_atom_site.auth_seq_id", "_atom_site.auth_comp_id",
        "_atom_site.auth_asym_id", "_atom_site.auth_atom_id", "_atom_site.pdbx_PDB_model_num"
    ]) + "\n" + "\n".join(rows) + "\n#\n"


def make_run(folder, ranks=(0, 1), with_request=True, names=("ALA", "CYS"), with_confidence=True):
    folder.mkdir()
    if with_request:
        (folder / "fold_synthetic_job_request.json").write_text(json.dumps(request()), encoding="utf-8")
    for rank in ranks:
        (folder / f"fold_synthetic_model_{rank}.cif").write_text(cif_text(names=names), encoding="utf-8")
        if with_confidence:
            (folder / f"fold_synthetic_summary_confidences_{rank}.json").write_text(json.dumps({"ptm": .8, "ranking_score": .9 - .1 * rank}), encoding="utf-8")
            (folder / f"fold_synthetic_full_data_{rank}.json").write_text(json.dumps({
                "atom_plddts": [20, 20, 80], "atom_chain_ids": ["A", "A", "A"],
                "pae": [[0, 2], [3, 0]], "token_chain_ids": ["A", "A"], "token_res_ids": [1, 2]
            }), encoding="utf-8")


def test_missing_download_is_waiting(tmp_path):
    result = import_alphafold(tmp_path / "absent", tmp_path / "retained", request())
    assert result["status"] == "WAITING_FOR_DOWNLOAD"
    assert result["models"] == []
    assert result["primary_model"] is None


def test_empty_and_readme_only_download_waits(tmp_path):
    source = tmp_path / "drop"
    source.mkdir()
    (source / "README.md").write_text("Drop manual downloads here.")
    assert ingest_files(source, tmp_path / "retained")["status"] == "WAITING_FOR_DOWNLOAD"


def test_nested_zip_gzip_preserves_originals_idempotently(tmp_path):
    source = tmp_path / "downloads"
    source.mkdir()
    compressed = gzip.compress(b">SYNTHETIC_QUERY\nAC\n>SYNTHETIC_HIT\nAC\n")
    nested = io.BytesIO()
    with zipfile.ZipFile(nested, "w") as archive:
        archive.writestr("alignment.a3m.gz", compressed)
    with zipfile.ZipFile(source / "manual.zip", "w") as archive:
        archive.writestr("nested/server.zip", nested.getvalue())
    first = ingest_files(source, tmp_path / "retained")
    second = ingest_files(source, tmp_path / "retained")
    assert first == second
    assert first["status"] == "IMPORTED"
    assert len(first["originals"]) == 1
    alignments = [item for item in first["files"] if item["path"].endswith(".a3m")]
    assert len(alignments) == 1
    assert len(alignments[0]["archive_chain"]) == 3
    assert (source / "manual.zip").read_bytes() == Path(first["originals"][0]["path"]).read_bytes()


@pytest.mark.parametrize("name", ["../escape.txt", "/absolute.txt", "C:/escape.txt", "folder\\..\\escape.txt", "file:stream", "CON.txt"])
def test_zip_traversal_and_special_paths_rejected(tmp_path, name):
    source = tmp_path / "manual.zip"
    with zipfile.ZipFile(source, "w") as archive:
        archive.writestr(name, "SYNTHETIC")
    result = ingest_files(source, tmp_path / "retained")
    assert result["status"] == "BLOCKED"
    assert len(result["originals"]) == 1
    assert not (tmp_path / "escape.txt").exists()


def test_zip_symlink_rejected(tmp_path):
    source = tmp_path / "manual.zip"
    with zipfile.ZipFile(source, "w") as archive:
        info = zipfile.ZipInfo("link")
        info.create_system = 3
        info.external_attr = (stat.S_IFLNK | 0o777) << 16
        archive.writestr(info, "../../elsewhere")
    result = ingest_files(source, tmp_path / "retained")
    assert result["status"] == "BLOCKED"
    assert "link/special" in result["issues"][0]


@pytest.mark.parametrize("tar_type", [tarfile.SYMTYPE, tarfile.LNKTYPE, tarfile.FIFOTYPE])
def test_tar_links_and_devices_rejected(tmp_path, tar_type):
    source = tmp_path / "manual.tar"
    with tarfile.open(source, "w") as archive:
        info = tarfile.TarInfo("link")
        info.type = tar_type
        info.linkname = "../outside"
        archive.addfile(info)
    result = ingest_files(source, tmp_path / "retained")
    assert result["status"] == "BLOCKED"


def test_tar_gzip_nested_accepted(tmp_path):
    source = tmp_path / "manual.tar.gz"
    with tarfile.open(source, "w:gz") as archive:
        info = tarfile.TarInfo("real_final.a3m")
        content = b">SYNTHETIC\nAC\n"
        info.size = len(content)
        archive.addfile(info, io.BytesIO(content))
    result = ingest_files(source, tmp_path / "retained")
    assert result["status"] == "IMPORTED"
    assert any(item["path"].endswith(".a3m") for item in result["files"])


def test_decompression_size_guard(tmp_path):
    source = tmp_path / "manual.gz"
    source.write_bytes(gzip.compress(b"x" * 10000))
    result = ingest_files(source, tmp_path / "retained", {"max_member_bytes": 200})
    assert result["status"] == "BLOCKED"
    assert len(result["originals"]) == 1
    assert not any(Path(item["path"]).name == "manual" for item in result["files"])


def test_archive_depth_guard(tmp_path):
    source = tmp_path / "manual.gz"
    source.write_bytes(gzip.compress(gzip.compress(b"SYNTHETIC")))
    result = ingest_files(source, tmp_path / "retained", {"max_depth": 1})
    assert result["status"] == "BLOCKED"
    assert "depth" in result["issues"][0]


def test_archive_count_guard(tmp_path):
    source = tmp_path / "manual.zip"
    with zipfile.ZipFile(source, "w") as archive:
        for index in range(5):
            archive.writestr(f"{index}.txt", "SYNTHETIC")
    assert ingest_files(source, tmp_path / "retained", {"max_members": 3})["status"] == "BLOCKED"


def test_duplicate_case_insensitive_archive_members_blocked(tmp_path):
    source = tmp_path / "manual.zip"
    with zipfile.ZipFile(source, "w") as archive:
        archive.writestr("same.txt", "SYNTHETIC_A")
        archive.writestr("SAME.txt", "SYNTHETIC_B")
    assert ingest_files(source, tmp_path / "retained")["status"] == "BLOCKED"


def test_tampered_retained_original_is_not_overwritten(tmp_path):
    source = tmp_path / "alignment.a3m"
    source.write_text(">SYNTHETIC\nAC\n")
    first = ingest_files(source, tmp_path / "retained")
    retained = Path(first["originals"][0]["path"])
    retained.write_text("TAMPERED")
    second = ingest_files(source, tmp_path / "retained")
    assert second["status"] == "BLOCKED"
    assert retained.read_text() == "TAMPERED"
    assert source.read_text() == ">SYNTHETIC\nAC\n"


def test_output_cannot_be_nested_under_source(tmp_path):
    with pytest.raises(ValueError, match="outside"):
        ingest_files(tmp_path, tmp_path / "retained")


def test_correct_request_and_seed_missing_or_changed():
    expected = request()
    assert verify_job_request(request(), expected)["status"] == "VERIFIED"
    missing = request()
    missing[0].pop("modelSeeds")
    assert verify_job_request(missing, expected)["status"] == "UNVERIFIED"
    assert verify_job_request(request(modelSeeds=["43"]), expected)["status"] == "BLOCKED"


@pytest.mark.parametrize("field,value", [("sequence", "CA"), ("count", 2), ("useStructureTemplate", True)])
def test_wrong_target_or_settings_blocks(field, value):
    returned = request()
    returned[0]["sequences"][0]["proteinChain"][field] = value
    assert verify_job_request(returned, request())["status"] == "BLOCKED"


def test_current_inventory_never_reuses_deleted_source(tmp_path):
    source = tmp_path / "drop"
    source.mkdir()
    old = source / "old.a3m"
    old.write_text(">SYNTHETIC_OLD\nAC\n")
    first = ingest_files(source, tmp_path / "retained")
    retained = Path(first["originals"][0]["path"])
    old.unlink()
    result = ingest_files(source, tmp_path / "retained")
    assert retained.exists()
    assert result["files"] == []
    assert result["status"] == "WAITING_FOR_DOWNLOAD"


def test_returned_embedded_auto_msa_is_preserved_but_provenance_unverified(tmp_path):
    source = tmp_path / "run1"
    make_run(source)
    returned = request()
    msa = ">SYNTHETIC_QUERY\nAC\n>SYNTHETIC_HIT\nAtC\n"
    returned[0]["sequences"][0]["proteinChain"]["unpairedMsa"] = msa
    (source / "fold_synthetic_job_request.json").write_text(json.dumps(returned))
    result = import_alphafold(source, tmp_path / "retained", request())
    assert result["status"] == "IMPORTED_UNVERIFIED"
    assert len(result["returned_alignments"]) == 1
    alignment = result["returned_alignments"][0]
    assert Path(alignment["path"]).read_text() == msa
    assert alignment["source_json_sha256"]


def test_missing_template_setting_stays_unverified():
    returned = request()
    returned[0]["sequences"][0]["proteinChain"].pop("useStructureTemplate")
    assert verify_job_request(returned, request())["status"] == "UNVERIFIED"


def test_missing_count_is_unverified_and_boolean_count_invalid():
    returned = request()
    returned[0]["sequences"][0]["proteinChain"].pop("count")
    assert verify_job_request(returned, request())["status"] == "UNVERIFIED"
    returned[0]["sequences"][0]["proteinChain"]["count"] = True
    assert verify_job_request(returned, request())["status"] == "BLOCKED"


def test_null_protein_entity_is_not_verified():
    returned = request()
    returned[0]["sequences"][0]["proteinChain"] = None
    assert verify_job_request(returned, request())["status"] == "BLOCKED"


def test_custom_msa_exact_evidence():
    expected = request()
    expected[0]["sequences"][0]["proteinChain"]["unpairedMsa"] = ">SYNTHETIC_QUERY\nAC\n>SYNTHETIC_HIT\nAtC\n"
    assert verify_job_request(expected, expected)["status"] == "VERIFIED"
    assert verify_job_request(request(), expected)["status"] == "UNVERIFIED"
    returned = request()
    returned[0]["sequences"][0]["proteinChain"]["unpairedMsa"] = ">SYNTHETIC_QUERY\nAC\n"
    assert verify_job_request(returned, expected)["status"] == "BLOCKED"
    assert verify_job_request(returned, request())["status"] == "UNVERIFIED"


def test_public_server_import_ranking_and_confidence(tmp_path):
    source = tmp_path / "run1"
    make_run(source, ranks=(0, 1, 2))
    result = import_alphafold(source, tmp_path / "retained", request())
    assert result["status"] == "VALIDATED"
    assert result["comparison_eligible"]
    assert len(result["models"]) == 3  # Never assume five models.
    assert result["primary_model"].endswith("_model_0.cif")
    confidence = result["models"][0]["confidence"]
    assert confidence["atom_mean_plddt"] == 40
    assert confidence["residue_mean_plddt"] == 50
    assert confidence["ptm"] == .8
    assert confidence["pae"] == [[0, 2], [3, 0]]
    assert confidence["pae_mapping_status"] == "VERIFIED"
    assert confidence["residue_plddt"][0]["auth_seq_id"] == "10"
    assert len(confidence["atom_mapping"]) == 3
    assert result["automatic_msa_status"].startswith("UNAVAILABLE")


def test_documented_filename_rank_when_scores_absent(tmp_path):
    source = tmp_path / "run1"
    make_run(source, ranks=(0, 1), with_confidence=False)
    result = import_alphafold(source, tmp_path / "retained", request())
    assert result["primary_model"].endswith("_model_0.cif")
    assert "filename" in result["primary_selection"] or "*_model_N.cif" in result["primary_selection"]
    assert result["models"][0]["confidence"]["residue_mean_plddt"] is None


def test_missing_rank_zero_requires_resolution_without_assuming_model_count(tmp_path):
    source = tmp_path / "run1"
    make_run(source, ranks=(1, 2))
    result = import_alphafold(source, tmp_path / "retained", request())
    assert result["primary_model"] is None
    assert len(result["models"]) == 2
    assert any("rank 0" in issue for issue in result["issues"])


def test_wrong_target_cif_blocks_even_when_request_matches(tmp_path):
    source = tmp_path / "run1"
    make_run(source, names=("ALA", "ALA"))
    result = import_alphafold(source, tmp_path / "retained", request())
    assert result["status"] == "BLOCKED"
    assert not result["comparison_eligible"]
    assert len(result["models"]) == 2


def test_missing_request_retains_models_unverified(tmp_path):
    source = tmp_path / "run1"
    make_run(source, with_request=False)
    result = import_alphafold(source, tmp_path / "retained", request())
    assert result["status"] == "IMPORTED_UNVERIFIED"
    assert result["verification"]["status"] == "UNVERIFIED"
    assert not result["comparison_eligible"]
    assert result["primary_model"]


def test_ranking_missing_requires_recorded_override(tmp_path):
    source = tmp_path / "run1"
    source.mkdir()
    (source / "mystery.cif").write_text(cif_text())
    blocked = import_alphafold(source, tmp_path / "retained", request())
    assert blocked["primary_model"] is None
    result = import_alphafold(source, tmp_path / "retained", request(), {"path": "mystery.cif", "rationale": "SYNTHETIC manual rank record for test"})
    assert result["primary_model"].endswith("mystery.cif")
    assert "SYNTHETIC manual" in result["primary_selection"]


def test_ranking_conflict_requires_override(tmp_path):
    source = tmp_path / "run1"
    make_run(source)
    (source / "fold_synthetic_summary_confidences_1.json").write_text(json.dumps({"ranking_score": 1.1}))
    result = import_alphafold(source, tmp_path / "retained", request())
    assert result["primary_model"] is None
    assert any("conflicts" in issue for issue in result["issues"])


def test_unmapped_atoms_never_become_residue_mean(tmp_path):
    model = tmp_path / "model.cif"
    model.write_text(cif_text())
    result = parse_confidence(model, {"ptm": .7}, {"atom_plddts": [10, 20, 30], "atom_chain_ids": ["B"] * 3})
    assert result["atom_mean_plddt"] == 20
    assert result["residue_mean_plddt"] is None
    assert result["ptm"] == .7


def test_nonfinite_or_out_of_range_confidence_stays_missing(tmp_path):
    model = tmp_path / "model.cif"
    model.write_text(cif_text())
    result = parse_confidence(model, {"ptm": 2}, {"atom_plddts": [float("nan"), 1000], "pae": [[-1]]})
    assert result["ptm"] is None
    assert result["atom_mean_plddt"] is None
    assert result["pae"] is None


@pytest.mark.parametrize("chains,residue_ids", [
    (["B", "B"], [1, 2]),
    (["A", "A"], [1, 1]),
    (["A", "A"], [1, 999]),
    (["A", "A"], [0, 1]),
    (["A", "A"], [True, 2]),
    (["A"], [1]),
])
def test_pae_mapping_rejects_wrong_duplicate_or_out_of_range_tokens(tmp_path, chains, residue_ids):
    model = tmp_path / "model.cif"
    model.write_text(cif_text())
    result = parse_confidence(model, full={"pae": [[0, 2], [3, 0]], "token_chain_ids": chains, "token_res_ids": residue_ids})
    assert result["pae"] == [[0, 2], [3, 0]]  # Data retained, identity never invented.
    assert result["pae_mapping_status"] == "UNVERIFIED"


def test_pae_mapping_preserves_label_and_author_residue_identity(tmp_path):
    model = tmp_path / "model.cif"
    model.write_text(cif_text())
    label = parse_confidence(model, full={"pae": [[0, 2], [3, 0]], "token_chain_ids": ["A", "A"], "token_res_ids": [1, 2]})
    assert label["pae_mapping_status"] == "VERIFIED"
    assert label["pae_mapping_convention"] == "label"
    assert label["pae_token_mapping"][0]["auth_seq_id"] == "10"
    assert label["pae_token_mapping"][1]["label_seq_id"] == "2"
    author = parse_confidence(model, full={"pae": [[0, 2], [3, 0]], "token_chain_ids": ["A", "A"], "token_res_ids": [10, 11]})
    assert author["pae_mapping_status"] == "VERIFIED"
    assert author["pae_mapping_convention"] == "author"
    assert author["pae_token_mapping"][0]["label_seq_id"] == "1"
