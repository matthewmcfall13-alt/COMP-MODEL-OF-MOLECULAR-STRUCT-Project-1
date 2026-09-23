"""Synthetic orchestration/package fixtures; no study inputs or outputs are used."""
from contextlib import ExitStack, redirect_stdout
import hashlib
import io
import json
from pathlib import Path
import runpy
import tempfile
import unittest
from unittest.mock import patch
import zipfile

import numpy as np

from project1.metrics import _write_ca_pdb
from project1.uploaded_analysis import _package, run_uploaded, uploads_available


class UploadedAnalysisIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="project1_SYNTHETIC_orchestration_")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.sequence = "ACDEF"
        self.coords = np.array([[0, 0, 0], [3, 0, 0], [3, 4, 0], [3, 4, 5], [1, 2, 6]], dtype=float)
        self.write_bytes("targets/target_sequences.fasta", b">SYNTHETIC_fixture\r\nACDEF\r\n")
        self.write_bytes("msa/synthetic/auto.a3m", b"#5\t1\r\n>query\r\nACDEF\r\n>duplicate\r\nACDEF\r\n")
        self.write_bytes("msa/synthetic/hmmer.a3m", b">query\r\nACDEF\r\n")
        for name in ("reference", "auto", "hmmer"):
            path = self.root / "structures" / (name + ".pdb")
            path.parent.mkdir(parents=True, exist_ok=True)
            _write_ca_pdb(path, list(range(1, 6)), self.sequence, self.coords)
            path.write_bytes(path.read_bytes().replace(b"\r\n", b"\n").replace(b"\n", b"\r\n"))
        self.write_bytes("targets/target_selection.md", b"SYNTHETIC fixture only.\r\n")
        self.write_bytes("analysis/analysis_summary.md", b"SYNTHETIC fixture; no scientific conclusions.\r\n")
        self.write_bytes("project1/synthetic_code_marker.py", b"# SYNTHETIC package marker only\r\n")
        self.config = {
            "repo_root": ".", "uploaded_analysis": {"enabled": True,
                "source_fasta": "targets/target_sequences.fasta", "targets": {
                    "SYNTHETIC": {"fasta_record_id": "SYNTHETIC_fixture",
                        "reference_path": "structures/reference.pdb", "reference_chain": "A",
                        "conditions": {
                            "Automatic_MMseqs2": {"alignment_path": "msa/synthetic/auto.a3m", "model_path": "structures/auto.pdb"},
                            "Custom_HMMER": {"alignment_path": "msa/synthetic/hmmer.a3m", "model_path": "structures/hmmer.pdb"},
                        }}}}}
        self.config_path = self.root / "config.json"
        self.save_config()

    def write_bytes(self, relative, content):
        path = self.root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
        return path

    def save_config(self):
        self.config_path.write_text(json.dumps(self.config), encoding="utf-8")

    def run_analysis(self):
        # Exercise real input parsing, all-row statistics, coordinate mapping,
        # scoring, reporting and integrity checks; renderer behavior is tested
        # separately and stubbed here to keep orchestration fixtures small.
        with ExitStack() as stack:
            stack.enter_context(patch("project1.uploaded_analysis._git", return_value="SYNTHETIC_revision"))
            stack.enter_context(patch("project1.uploaded_analysis._plot_alignment_comparison"))
            stack.enter_context(patch("project1.uploaded_msa._figures", return_value={}))
            stack.enter_context(patch("project1.uploaded_structures._plots", return_value="SYNTHETIC_renderer_stub"))
            stack.enter_context(redirect_stdout(io.StringIO()))
            return run_uploaded(self.config_path, package=False)

    def test_complete_analysis_records_inventory_without_mutating_sources(self):
        before = {p: p.read_bytes() for p in self.root.rglob("*") if p.is_file()}
        result = self.run_analysis()
        self.assertEqual(result["status"], "COMPLETED_FOR_AVAILABLE_UPLOADS")
        self.assertEqual(result["originals_integrity"], "UNCHANGED")
        self.assertEqual(result["changed_inputs"], [])
        self.assertTrue(result["controls_status"].startswith("UNVERIFIED"))
        self.assertEqual([row["status"] for row in result["alignments"]], ["NO_NON_QUERY_DIVERSITY", "QUERY_ONLY"])
        self.assertEqual(len(result["structures"]["rows"]), 2)
        for row in result["structures"]["rows"]:
            self.assertLess(row["rmsd_ca_angstrom"], 1e-12)
            self.assertEqual(row["lddt_ca_common_mask"], 1)
            self.assertEqual(row["experimental_controls_status"], "UNVERIFIED")
        for item in result["input_inventory"]:
            if item["exists"]:
                self.assertEqual(item["sha256"], hashlib.sha256(before[Path(item["path"])]).hexdigest())
        self.assertEqual(before, {p: p.read_bytes() for p in before})
        output = Path(result["output_directory"])
        self.assertTrue((output / "REPORT.md").is_file())
        pointer = json.loads((self.root / "outputs/uploaded_analysis_current.json").read_text())
        self.assertEqual(pointer["output_directory"], str(output))
        self.assertEqual(json.loads((output / "results.json").read_text())["input_inventory"], result["input_inventory"])

    def test_missing_alignment_marks_partial_and_keeps_other_condition(self):
        self.config["uploaded_analysis"]["targets"]["SYNTHETIC"]["conditions"]["Custom_HMMER"]["alignment_path"] = "msa/synthetic/never_created.a3m"
        self.save_config()
        result = self.run_analysis()
        self.assertEqual(result["status"], "PARTIAL_REQUIRES_REVIEW")
        self.assertEqual(result["originals_integrity"], "UNCHANGED")
        self.assertEqual([row["condition"] for row in result["alignments"]], ["Automatic_MMseqs2"])
        self.assertTrue(any("Custom_HMMER" in issue and "FileNotFoundError" in issue for issue in result["issues"]))
        missing = [item for item in result["input_inventory"] if item["relative_path"].endswith("never_created.a3m")]
        self.assertEqual(len(missing), 1)
        self.assertFalse(missing[0]["exists"])
        self.assertIsNone(missing[0]["sha256"])
        self.assertTrue((Path(result["output_directory"]) / "diagnostics/SYNTHETIC_Custom_HMMER.txt").is_file())

    def test_package_preserves_scientific_bytes_and_strips_notebook_outputs(self):
        result = self.run_analysis()
        notebook = {"nbformat": 4, "nbformat_minor": 5, "metadata": {}, "cells": [{
            "cell_type": "code", "metadata": {}, "source": ["print('SYNTHETIC')"], "execution_count": 9,
            "outputs": [{"output_type": "stream", "name": "stdout", "text": "SYNTHETIC_TRANSIENT_OUTPUT"}]}]}
        notebook_path = self.write_bytes("Project1_Server_Workflow.ipynb", json.dumps(notebook).encode())
        before = notebook_path.read_bytes()
        archive_path = _package(self.root, Path(result["output_directory"]), result["input_inventory"])
        with zipfile.ZipFile(archive_path) as archive:
            self.assertIsNone(archive.testzip())
            for item in result["input_inventory"]:
                if item["exists"] and Path(item["path"]).suffix.lower() in {".a3m", ".fasta", ".pdb"}:
                    content = archive.read(item["relative_path"])
                    self.assertEqual(content, Path(item["path"]).read_bytes())
                    self.assertEqual(hashlib.sha256(content).hexdigest(), item["sha256"])
            packaged_notebook = json.loads(archive.read("Project1_Server_Workflow.ipynb"))
            self.assertEqual(packaged_notebook["cells"][0]["outputs"], [])
            self.assertIsNone(packaged_notebook["cells"][0]["execution_count"])
            self.assertEqual(packaged_notebook["cells"][0]["source"], notebook["cells"][0]["source"])
        self.assertEqual(notebook_path.read_bytes(), before)

    def test_auto_detection_requires_enabled_configuration_and_fasta(self):
        self.assertTrue(uploads_available(self.config_path))
        self.config["uploaded_analysis"]["enabled"] = False
        self.save_config()
        self.assertFalse(uploads_available(self.config_path))
        self.config["uploaded_analysis"]["enabled"] = True
        self.config["uploaded_analysis"]["source_fasta"] = "targets/never_created.fasta"
        self.save_config()
        self.assertFalse(uploads_available(self.config_path))

    def test_cli_routes_auto_and_explicit_modes_without_running_jobs(self):
        scenarios = [
            (["--mode", "auto"], True, "uploaded"),
            (["--mode", "auto"], False, "manual"),
            (["--mode", "uploaded"], False, "uploaded"),
            (["--mode", "manual-server"], True, "manual"),
        ]
        for options, available, expected in scenarios:
            with self.subTest(options=options, available=available), ExitStack() as stack:
                stack.enter_context(patch("sys.argv", ["project1", "--config", str(self.config_path), "--no-package", *options]))
                stack.enter_context(patch("project1.uploaded_analysis.uploads_available", return_value=available))
                uploaded = stack.enter_context(patch("project1.uploaded_analysis.run_uploaded"))
                manual = stack.enter_context(patch("project1.workflow.run"))
                stack.enter_context(patch("project1.workflow.prepare"))
                runpy.run_module("project1.__main__", run_name="__main__")
                chosen, other = (uploaded, manual) if expected == "uploaded" else (manual, uploaded)
                chosen.assert_called_once_with(str(self.config_path), package=False)
                other.assert_not_called()

    def test_auto_missing_current_fasta_does_not_prepare_legacy_study(self):
        self.config["uploaded_analysis"]["source_fasta"] = "targets/never_created.fasta"
        self.save_config()
        with patch("sys.argv", ["project1", "--config", str(self.config_path)]), patch("project1.workflow.run") as manual:
            with self.assertRaisesRegex(FileNotFoundError, "Configured AF2 study FASTA is missing"):
                runpy.run_module("project1.__main__", run_name="__main__")
            manual.assert_not_called()
        self.assertFalse((self.root / "START_HERE.md").exists())


if __name__ == "__main__":
    unittest.main()
