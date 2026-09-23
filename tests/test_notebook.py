"""Fresh-kernel FASTA-only smoke test in an isolated copy, never real output folders."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys

import nbformat
from nbclient import NotebookClient
from jupyter_client import KernelManager


def test_fresh_kernel_notebook_fasta_only(tmp_path):
    root = Path(__file__).resolve().parents[1]
    workspace = tmp_path / "isolated_notebook"
    workspace.mkdir()
    for name in ("project1", "docs"):
        shutil.copytree(root / name, workspace / name, ignore=shutil.ignore_patterns("__pycache__"))
    for name in ("config.json", "Project1_Server_Workflow.ipynb", "T1106s1.fasta", "T1122.fasta", "T1151s2.fasta", "label.xlsx"):
        shutil.copy2(root / name, workspace / name)
    # The completed AF2 study uses Nancy's targets in the main configuration.
    # Exercise the preserved initial FASTA-only workflow with its explicit
    # archived configuration, without relabeling old inputs as current results.
    legacy = json.loads((root / "docs/legacy_af3_config.json").read_text())
    legacy["repo_root"] = "."
    (workspace / "config.json").write_text(json.dumps(legacy), encoding="utf-8")
    before = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in workspace.glob("*.fasta")}
    nb = nbformat.read(workspace / "Project1_Server_Workflow.ipynb", as_version=4)
    assert nb.cells and all(c.cell_type == "code" for c in nb.cells)
    km = KernelManager(kernel_name="python3")
    km.kernel_spec.argv[0] = sys.executable
    env = dict(os.environ, PROJECT1_CONFIG=str(workspace / "config.json"))
    NotebookClient(nb, km=km, timeout=120, resources={"metadata": {"path": str(workspace)}}).execute(env=env)
    assert not any(output.output_type == "error" for cell in nb.cells for output in cell.outputs)
    status = json.loads((workspace / "outputs/status.json").read_text())
    assert len(status["targets"]) == 3
    assert all(row["AUTO"]["state"] == "PREPARED" for row in status["targets"].values())
    assert all(row["DEEPMSA"]["state"] == "WAITING_FOR_DEEPMSA" for row in status["targets"].values())
    assert before == {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in workspace.glob("*.fasta")}
    assert len(list((workspace / "outputs/reviews").glob("*.zip"))) == 1
    assert len((workspace / "outputs/accuracy_summary.csv").read_text().splitlines()) == 1
    assert len((workspace / "outputs/model_confidence.csv").read_text().splitlines()) == 1
