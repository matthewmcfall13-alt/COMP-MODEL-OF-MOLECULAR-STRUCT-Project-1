"""Synthetic coordinate audit fixtures, generated only in temporary folders."""
import csv
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from project1.metrics import _write_ca_pdb
from project1.uploaded_structures import audit_structures, rounding_comparison


class UploadedStructureAuditTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="project1_SYNTHETIC_uploaded_")
        self.root = Path(self.temporary.name)
        self.sequence = "ACDEF"
        self.coords = np.array([[0, 0, 0], [3, 0, 0], [3, 4, 0], [3, 4, 5], [1, 2, 6]], dtype=float)
        self.write_model("reference.pdb", range(1, 6), self.coords)
        self.write_model("auto.pdb", range(1, 6), self.coords)
        rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        self.write_model("custom.pdb", range(1, 6), self.coords @ rotation + [9, -12, 7])
        self.targets = {"SYNTHETIC": {"sequence": self.sequence, "reference_path": "reference.pdb",
            "conditions": {"AUTO": {"model_path": "auto.pdb"}, "HMMER": {"model_path": "custom.pdb"}}}}

    def tearDown(self):
        self.temporary.cleanup()

    def write_model(self, filename, positions, coordinates, bfactor=80):
        path = self.root / filename
        _write_ca_pdb(path, list(positions), self.sequence, coordinates)
        path.write_text("\n".join(line[:60] + f"{bfactor:6.2f}" + line[66:] if line.startswith("ATOM") else line
                                  for line in path.read_text().splitlines()) + "\n", encoding="ascii")

    def audit(self):
        # Rendering is separately exercised by the real-data notebook; these
        # synthetic unit tests focus on coordinate/mapping correctness.
        with patch("project1.uploaded_structures._plots", return_value=None):
            return audit_structures(self.root, self.root / "derived", self.targets)

    def test_self_rigid_fit_profiles_preserve_sources_and_json(self):
        before = {p: p.read_bytes() for p in self.root.glob("*.pdb")}
        result = self.audit()
        self.assertEqual(result["status"], "COMPLETE")
        self.assertEqual(len(result["rows"]), 2)
        for row in result["rows"]:
            self.assertLess(row["rmsd_ca_angstrom"], 1e-12)
            self.assertEqual(row["lddt_ca_common_mask"], 1)
            self.assertEqual(row["ca_bfactor_residue_mean_all_mapped"], 80)
            self.assertEqual(row["experimental_controls_status"], "UNVERIFIED")
            self.assertEqual(row["coordinate_comparison_status"], "COMPARABLE")
            self.assertIsNone(row["reported_ptm"])
        self.assertEqual(before, {p: p.read_bytes() for p in before})
        self.assertEqual(json.loads((self.root / "derived/result.json").read_text()), result)
        with (self.root / "derived/SYNTHETIC/per_residue_profiles.csv").open() as handle:
            profiles = list(csv.DictReader(handle))
        self.assertEqual(len(profiles), 10)
        self.assertEqual(profiles[0]["reference_residue_id"], "1")
        self.assertEqual(len(result["targets"]["SYNTHETIC"]["inputs"]["AUTO"]["sha256"]), 64)

    def test_frozen_mask_uses_every_condition_and_not_distance(self):
        self.write_model("third.pdb", [1, 2, 3, 4], self.coords[:4])
        self.targets["SYNTHETIC"]["conditions"]["DEEPMSA"] = {"model_path": "third.pdb"}
        outlier = self.coords.copy()
        outlier[3] += [100, 30, 2]
        self.write_model("auto.pdb", range(1, 6), outlier)
        result = self.audit()
        self.assertEqual(result["targets"]["SYNTHETIC"]["frozen_mask_positions"], [1, 2, 3, 4])
        self.assertTrue(all(row["scored_residue_count"] == 4 for row in result["rows"]))
        self.assertGreater(result["rows"][0]["rmsd_ca_angstrom"], 10)
        self.assertLess(result["rows"][0]["lddt_ca_common_mask"], 1)
        self.assertEqual(result["rows"][0]["scored_fraction_of_target"], 0.8)

    def test_ambiguous_mapping_blocks_entire_target(self):
        self.sequence = "AAAAA"
        self.targets["SYNTHETIC"]["sequence"] = self.sequence
        for filename in ("reference.pdb", "auto.pdb"):
            self.write_model(filename, range(1, 6), self.coords)
        self.write_model("custom.pdb", [1, 2, 3, 4], self.coords[:4])
        result = self.audit()
        self.assertEqual(result["status"], "PARTIAL")
        self.assertEqual(result["rows"], [])
        self.assertIn("Ambiguous", result["targets"]["SYNTHETIC"]["diagnostic"])

    def test_explicit_identity_checked_mapping_resolves_ambiguity(self):
        self.sequence = "AAAAA"
        self.targets["SYNTHETIC"]["sequence"] = self.sequence
        for filename in ("reference.pdb", "auto.pdb"):
            self.write_model(filename, range(1, 6), self.coords)
        self.write_model("custom.pdb", [1, 2, 3, 4], self.coords[:4])
        self.targets["SYNTHETIC"]["conditions"]["HMMER"]["mapping"] = {i: str(i) for i in range(1, 5)}
        result = self.audit()
        self.assertEqual(result["status"], "COMPLETE")
        self.assertEqual(result["rows"][0]["scored_residue_count"], 4)

    def test_wrong_chain_or_too_few_residues_block_all(self):
        self.targets["SYNTHETIC"]["conditions"]["HMMER"]["model_chain"] = "X"
        self.assertEqual(self.audit()["rows"], [])
        self.targets["SYNTHETIC"]["conditions"]["HMMER"].pop("model_chain")
        self.write_model("custom.pdb", [1, 2], self.coords[:2])
        result = self.audit()
        self.assertEqual(result["rows"], [])
        self.assertIn("Fewer than three", result["targets"]["SYNTHETIC"]["diagnostic"])

    def test_reported_metrics_are_compared_but_tm_and_ptm_unverified(self):
        (self.root / "analysis").mkdir()
        (self.root / "analysis/casp15_structure_comparison_metrics.csv").write_text(
            "Target,Condition,RMSD_CA_A,lDDT_CA,Mapped_residues,pLDDT,pTM,TM_score\n"
            "SYNTHETIC,Automatic MMseqs2,0.00,100.00,5,80.00,0.99,0.98\n"
            "SYNTHETIC,Custom HMMER,1.00,99.00,5,85.00,0.98,0.97\n", encoding="utf-8")
        result = self.audit()
        auto, hmmer = result["rows"]
        self.assertEqual(auto["rmsd_report_comparison"], "MATCH_WITHIN_ROUNDING")
        self.assertEqual(auto["lddt_report_comparison"], "MATCH_WITHIN_ROUNDING")
        self.assertEqual(auto["reported_ptm"], "0.99")
        self.assertTrue(auto["ptm_verification"].startswith("UNVERIFIED"))
        self.assertTrue(auto["tm_score_verification"].startswith("UNVERIFIED"))
        self.assertEqual(hmmer["rmsd_report_comparison"], "DIFFERS")
        self.assertEqual(hmmer["lddt_report_comparison"], "DIFFERS")
        self.assertEqual(hmmer["plddt_report_comparison"], "DIFFERS_AGGREGATION_UNVERIFIED")

    def test_rounding_tolerances_and_invalid_values(self):
        self.assertEqual(rounding_comparison(0.8749, "0.87")["status"], "MATCH_WITHIN_ROUNDING")
        self.assertEqual(rounding_comparison(0.8751, "0.87")["status"], "DIFFERS")
        self.assertEqual(rounding_comparison(96.884, "96.88")["tolerance"], 0.005)
        self.assertEqual(rounding_comparison(None, "96.88")["status"], "UNAVAILABLE")
        self.assertEqual(rounding_comparison(1, "NaN")["status"], "INVALID_REPORTED_VALUE")

    def test_reported_mask_count_mismatch_not_comparable(self):
        (self.root / "analysis").mkdir()
        (self.root / "analysis/casp15_structure_comparison_metrics.csv").write_text(
            "Target,Condition,RMSD_CA_A,lDDT_CA,Mapped_residues\n"
            "SYNTHETIC,Automatic MMseqs2,0.00,100.00,4\n", encoding="utf-8")
        row = self.audit()["rows"][0]
        self.assertEqual(row["rmsd_report_comparison"], "NOT_COMPARABLE_MASK_COUNT")
        self.assertEqual(row["lddt_report_comparison"], "NOT_COMPARABLE_MASK_COUNT")

    def _declare_known_reference(self, pdb_id, experimental_method=None):
        self.targets["T1188"] = self.targets.pop("SYNTHETIC")
        path = self.root / "reference.pdb"
        header = "HEADER".ljust(62) + pdb_id.ljust(4) + "\n"
        if experimental_method is not None:
            header += "EXPDTA    " + experimental_method + "\n"
        path.write_text(header + path.read_text(encoding="ascii"), encoding="ascii")

    def test_known_assignment_rejects_wrong_header_id_despite_matching_sequence(self):
        self._declare_known_reference("8BBT", "X-RAY DIFFRACTION")
        result = self.audit()
        self.assertEqual(result["rows"], [])
        target = result["targets"]["T1188"]
        self.assertEqual(target["status"], "BLOCKED")
        self.assertIn("does not match documented CASP15 assignment 8C6Z", target["diagnostic"])
        self.assertFalse(target["casp_assignment"]["uploaded_header_matches_assignment"])

    def test_known_assignment_requires_xray_experimental_evidence(self):
        self._declare_known_reference("8C6Z", "THEORETICAL MODEL")
        result = self.audit()
        self.assertEqual(result["rows"], [])
        self.assertIn("lacks the required experimental", result["targets"]["T1188"]["diagnostic"])
        reference = self.root / "reference.pdb"
        reference.write_text("\n".join(line for line in reference.read_text().splitlines()
                                       if not line.startswith("EXPDTA")) + "\n", encoding="ascii")
        missing_evidence = self.audit()
        self.assertEqual(missing_evidence["rows"], [])
        self.assertEqual(missing_evidence["targets"]["T1188"]["status"], "BLOCKED")

    def test_known_assignment_with_matching_header_and_xray_evidence_scores(self):
        self._declare_known_reference("8C6Z", "X-RAY DIFFRACTION")
        result = self.audit()
        self.assertEqual(result["status"], "COMPLETE")
        self.assertEqual(len(result["rows"]), 2)
        self.assertTrue(result["targets"]["T1188"]["casp_assignment"]["uploaded_header_matches_assignment"])

    def test_invalid_confidence_not_interpreted_as_plddt(self):
        self.write_model("auto.pdb", range(1, 6), self.coords, bfactor=110)
        row = self.audit()["rows"][0]
        self.assertFalse(row["ca_bfactor_0_100_valid"])
        self.assertIsNone(row["ca_bfactor_residue_mean_all_mapped"])


if __name__ == "__main__":
    unittest.main()
