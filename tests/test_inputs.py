"""SYNTHETIC fixtures; none are project targets or experimental results."""
import json
from pathlib import Path

import pytest

from project1.inputs import discover_fasta, parse_a3m, read_fasta, select_a3m, server_request, verify_pair
from project1.workflow import prepare


def test_fasta_normalization_preserves_bytes(tmp_path):
    source = tmp_path / "original.fasta"
    original = b">synthetic demonstration\r\nac de\r\nFG\r\n"
    source.write_bytes(original)
    rec = read_fasta(source)[0]
    assert rec["sequence"] == "ACDEFG"
    assert not rec["errors"]
    assert rec["normalization_changed"]
    assert source.read_bytes() == original


@pytest.mark.parametrize("sequence", ["", "ACX", "AC-", "AC*", "ACU", "ACB", "ACZ", "ACJ", "ACO", "AC1"])
def test_invalid_protein_is_not_repaired(tmp_path, sequence):
    p = tmp_path / "bad.fa"
    p.write_text(">synthetic\n" + sequence + "\n")
    record = read_fasta(p)[0]
    assert record["errors"]
    assert record["sequence"] == sequence


def test_nonfasta_and_no_records(tmp_path):
    for text in ("ACDE", ""):
        p = tmp_path / "bad.fa"
        p.write_text(text)
        with pytest.raises(ValueError):
            read_fasta(p)


def test_multiple_records_duplicates_and_exclusion(tmp_path):
    (tmp_path / "input.fa").write_text(">a\nACD\n>b\nACE\n>a\nACF\n>c\nACE\n")
    (tmp_path / "tests").mkdir()
    (tmp_path / "tests" / "excluded.fa").write_text(">test_only\nAAA\n")
    records, _ = discover_fasta(tmp_path)
    assert len(records) == 4
    assert all(r["validation_status"] == "BLOCKED" for r in records)


def test_server_schema_and_a3m_exact_bytes():
    text = ">synthetic\nACDE\n>homolog\nACdD-\n"
    auto = server_request("synthetic", "ACDE")
    deep = server_request("synthetic", "ACDE", msa=text)
    assert isinstance(auto, list) and len(auto) == 1
    assert auto[0]["dialect"] == "alphafoldserver" and auto[0]["version"] == 1
    assert auto[0]["modelSeeds"] == ["42"]
    assert auto[0]["sequences"] == [{"proteinChain": {"sequence": "ACDE", "count": 1, "useStructureTemplate": False}}]
    assert deep[0]["sequences"][0]["proteinChain"]["unpairedMsa"] == text
    assert parse_a3m(text, "ACDE")["analysis_rows"] == ["ACDE", "ACD-"]
    assert verify_pair(auto, deep)
    deep[0]["modelSeeds"] = ["43"]
    with pytest.raises(ValueError, match="differ"):
        verify_pair(auto, deep)


@pytest.mark.parametrize("seed", [-1, 2**32, "42", 42.0, True])
def test_seed_validation(seed):
    with pytest.raises(ValueError):
        server_request("synthetic", "ACD", seed)


@pytest.mark.parametrize("text, error", [
    (">q\nACDF\n>h\nACDE\n", "first sequence"),
    (">q\nACDE\n>h\nACD\n", "widths"),
    (">q\nACDE\n", "Query-only"),
    (">q\nACdDE\n>h\nACDE\n", "first sequence"),
    (">q\nACDE\n>h\nAC*E\n", "unsupported"),
])
def test_a3m_blocks(text, error):
    with pytest.raises(ValueError, match=error):
        parse_a3m(text, "ACDE")


def test_a3m_ambiguous_candidates(tmp_path):
    files = [tmp_path / "small.a3m", tmp_path / "large.a3m"]
    with pytest.raises(ValueError, match="Multiple"):
        select_a3m(files)
    assert select_a3m(files, "small.a3m")[0] == files[0]
    with pytest.raises(ValueError, match="evidence"):
        select_a3m(files, service_final={"path": "small.a3m"})


def make_config(tmp_path):
    (tmp_path / "input.fasta").write_text(">synthetic\nACDE\n")
    config = {"repo_root": ".", "data_root": "data", "selected_targets": ["synthetic"], "pilot": "synthetic", "seed": 42}
    path = tmp_path / "config.json"
    path.write_text(json.dumps(config))
    return path


def test_prepare_idempotent_preserves_manual_and_waiting(tmp_path):
    path = make_config(tmp_path)
    _, first = prepare(path)
    manifest_path = tmp_path / "target_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["targets"]["synthetic"].update({"notes": "manual annotation", "reference_chain": "Z"})
    manifest_path.write_text(json.dumps(manifest))
    auto = Path(first["targets"]["synthetic"]["AUTO"]["request_path"])
    original = auto.read_bytes()
    _, second = prepare(path)
    assert auto.read_bytes() == original
    assert len(list((tmp_path / "submissions").rglob("*.json"))) == 1
    assert second["targets"]["synthetic"]["DEEPMSA"]["state"] == "WAITING_FOR_DEEPMSA"
    assert json.loads(manifest_path.read_text())["targets"]["synthetic"]["notes"] == "manual annotation"
    assert json.loads(manifest_path.read_text())["targets"]["synthetic"]["reference_chain"] == "Z"


def test_empty_selection_no_arbitrary_pick(tmp_path):
    path = make_config(tmp_path)
    config = json.loads(path.read_text()); config["selected_targets"] = []
    path.write_text(json.dumps(config))
    _, state = prepare(path)
    assert not state["targets"]
    assert not (tmp_path / "submissions").exists()


def test_invalid_target_independent_valid_proceeds(tmp_path):
    path = make_config(tmp_path)
    (tmp_path / "invalid.fa").write_text(">bad\nAC*\n")
    config = json.loads(path.read_text()); config["selected_targets"].append("bad")
    path.write_text(json.dumps(config))
    _, state = prepare(path)
    assert state["targets"]["synthetic"]["AUTO"]["state"] == "PREPARED"
    assert state["targets"]["bad"]["AUTO"]["state"] == "BLOCKED"
    assert not (tmp_path / "submissions" / "bad").exists()
