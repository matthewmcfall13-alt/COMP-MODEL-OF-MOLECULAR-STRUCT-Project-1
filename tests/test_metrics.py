"""Synthetic unit fixtures only; never write into study/result folders."""
import json
import os
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from Bio.PDB import MMCIFIO, PDBParser

from project1.metrics import (
    MappingError, _write_ca_pdb, alignment_statistics, analyze_pair,
    kabsch_rmsd, lddt_ca, map_to_target, parse_tmscore_output,
    read_ca_structure, run_tmscore,
)


# This text is SYNTHETIC, not captured experimental data or a real tool run.
# Format/normalization contract checked against official TMscore.h output code.
SYNTHETIC_TMSCORE_OUTPUT = """
Structure1: model.pdb    Length=    5
Structure2: native.pdb   Length=    5 (by which all scores are normalized)
Number of residues in common=    5
RMSD of  the common residues=    0.125
TM-score    = 0.9876  (d0= 0.50)
GDT-TS-score= 0.8000 %(d<1)=0.6000 %(d<2)=0.8000 %(d<4)=0.8000 %(d<8)=1.0000
TM-score    = 0.4567  (if normalized by user-specified LN=10.00 and d0=0.50)
"""


class AlignmentStatisticsTests(unittest.TestCase):
    def test_column_counts_and_eligible_denominator(self):
        result = alignment_statistics(["AC", "AC", "A-"])
        # A/A=4, C/C=9. Three A pairs and one C pair: 21/4.
        self.assertEqual(result["blosum62_sum_of_pairs"], 21)
        self.assertEqual(result["blosum62_eligible_pair_denominator"], 4)
        self.assertEqual(result["blosum62_normalized_sum_of_pairs"], 5.25)
        self.assertEqual(result["duplicate_rows_beyond_first"], 1)
        self.assertEqual(result["sequence_count_including_query"], 3)
        self.assertAlmostEqual(result["gap_fraction"], 1 / 6)
        self.assertEqual(result["positions"][1]["non_gap_coverage"], 2 / 3)
        self.assertEqual(result["mean_identity_to_query_excluding_query"], 1)

    def test_no_eligible_pairs_is_missing(self):
        result = alignment_statistics(["AC", "--"])
        self.assertIsNone(result["blosum62_normalized_sum_of_pairs"])
        self.assertIsNone(result["rows"][1]["identity_to_query"])
        self.assertEqual(result["neff"], 1)

    def test_noncanonical_homologs_excluded_from_blosum_not_coverage(self):
        result = alignment_statistics(["AC", "AX"])
        self.assertEqual(result["blosum62_eligible_pair_denominator"], 1)
        self.assertEqual(result["positions"][1]["non_gap_coverage"], 1)
        self.assertEqual(result["positions"][1]["canonical_coverage"], 0.5)

    def test_neff_exact_and_bounded(self):
        self.assertEqual(alignment_statistics(["AC", "AC"])["neff"], 1)
        result = alignment_statistics(["AC", "AC"], neff_max_rows=1)
        self.assertIsNone(result["neff"])
        self.assertEqual(result["neff_status"], "UNAVAILABLE_EXACT_ROW_LIMIT")
        self.assertFalse(result["neff_approximation"])

    def test_lowercase_insertions_must_be_removed_upstream(self):
        with self.assertRaisesRegex(ValueError, "lowercase"):
            alignment_statistics(["AC", "Ac"])
        with self.assertRaises(ValueError):
            alignment_statistics(["AC", "A"])


class CoordinateTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="project1_synthetic_metrics_")
        self.root = Path(self.tmp.name)
        self.sequence = "ACDEF"
        self.coords = np.array([[0, 0, 0], [3, 0, 0], [3, 4, 0], [3, 4, 5], [1, 2, 6]], dtype=float)
        self.ref = self.root / "reference_SYNTHETIC.pdb"
        _write_ca_pdb(self.ref, list(range(1, 6)), self.sequence, self.coords)

    def tearDown(self):
        self.tmp.cleanup()

    def test_self_and_rigid_transform_rmsd(self):
        self.assertLess(kabsch_rmsd(self.coords, self.coords)["rmsd_angstrom"], 1e-12)
        angle = 0.7
        rotation = np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
        moved = self.coords @ rotation + [8, -9, 22]
        self.assertLess(kabsch_rmsd(moved, self.coords)["rmsd_angstrom"], 1e-12)
        self.assertEqual(lddt_ca(moved, self.coords)["score"], 1)
        perturbed = moved.copy()
        perturbed[-1] += [8, -4, 2]
        self.assertGreater(kabsch_rmsd(perturbed, self.coords)["rmsd_angstrom"], 1)
        self.assertLess(lddt_ca(perturbed, self.coords)["score"], 1)

    def test_lddt_exact_thresholds_and_missing_mobile(self):
        reference = np.array([[0, 0, 0], [3, 0, 0]], dtype=float)
        mobile = np.array([[0, 0, 0], [3.5, 0, 0]], dtype=float)
        result = lddt_ca(mobile, reference)
        self.assertEqual(result["conserved_by_threshold"], [0, 1, 1, 1])
        self.assertEqual(result["score"], 0.75)
        mobile[-1] = np.nan
        self.assertEqual(lddt_ca(mobile, reference)["score"], 0)
        self.assertIsNone(lddt_ca(reference * 10, reference * 10)["score"])

    def test_no_reflection_fit(self):
        mirrored = self.coords.copy()
        mirrored[:, 0] *= -1
        fit = kabsch_rmsd(mirrored, self.coords)
        self.assertGreater(fit["rmsd_angstrom"], 0.1)
        self.assertAlmostEqual(np.linalg.det(fit["rotation_row_vector"]), 1)

    def test_pdb_cif_chain_number_and_insertion_mapping(self):
        original = self.ref.read_text()
        modified = []
        for index, line in enumerate(original.splitlines()):
            if line.startswith("ATOM"):
                number = 41 + index
                insertion = "A" if index == 1 else " "
                line = line[:21] + "B" + f"{number:4d}" + insertion + line[27:]
            modified.append(line)
        self.ref.write_text("\n".join(modified) + "\n")
        structure = read_ca_structure(self.ref, "B")
        self.assertEqual(structure["residues"][1]["residue_id"], "42A")
        mapping = map_to_target(self.sequence, structure)
        self.assertEqual(mapping[2]["number"], 42)
        with self.assertRaises(MappingError):
            read_ca_structure(self.ref, "A")
        cif = self.root / "reference_SYNTHETIC.cif"
        writer = MMCIFIO()
        writer.set_structure(PDBParser(QUIET=True).get_structure("synthetic", self.ref))
        writer.save(str(cif))
        cif_structure = read_ca_structure(cif, "B")
        self.assertEqual(cif_structure["residues"][1]["insertion_code"], "A")
        np.testing.assert_allclose(cif_structure["residues"][1]["xyz"], self.coords[1])

    def test_ambiguous_mapping_requires_explicit_resolution(self):
        path = self.root / "ambiguous_SYNTHETIC.pdb"
        _write_ca_pdb(path, [1, 2, 3], "AAA", self.coords[:3])
        structure = read_ca_structure(path, "A")
        with self.assertRaisesRegex(MappingError, "Ambiguous"):
            map_to_target("AAAA", structure)
        result = map_to_target("AAAA", structure, {1: "1", 3: "2", 4: "3"})
        self.assertEqual(sorted(result), [1, 3, 4])
        with self.assertRaisesRegex(MappingError, "mismatch"):
            map_to_target("ACAA", structure, {2: "2"})

    def test_missing_renumbered_residues_and_frozen_pair_mask(self):
        auto = self.root / "AUTO_SYNTHETIC.pdb"
        deep = self.root / "DEEPMSA_SYNTHETIC.pdb"
        _write_ca_pdb(auto, [1, 2, 3, 4, 5], self.sequence, self.coords + 10)
        _write_ca_pdb(deep, [1, 2, 4, 5], self.sequence, self.coords[[0, 1, 3, 4]] - 3)
        output = self.root / "synthetic_analysis"
        with patch("project1.metrics.shutil.which", return_value=None):
            result = analyze_pair(self.sequence, self.ref, "A", auto, deep, output, reference_length=5)
        self.assertEqual(result["status"], "ANALYZED", result)
        self.assertEqual(result["scored_residue_count"], 4)
        self.assertEqual(result["reference_normalization_length"], 5)
        self.assertAlmostEqual(result["scored_fraction_of_target"], 0.8)
        for condition in ("AUTO", "DEEPMSA"):
            self.assertLess(result["metrics"][condition]["rmsd_angstrom"], 1e-6)
            self.assertEqual(result["metrics"][condition]["lddt_ca_common_mask"], 1)
            self.assertIsNone(result["metrics"][condition]["tm_score"])
        mask = (output / "frozen_mask.csv").read_text()
        self.assertEqual(mask.splitlines()[1:], ["1", "2", "4", "5"])
        self.assertIn("matchNumbering true", (output / "overlay.cxc").read_text())
        self.assertNotIn("cutoffDistance ", "\n".join(line for line in (output / "overlay.cxc").read_text().splitlines() if not line.startswith("#")))
        self.assertEqual(json.loads((output / "metric_summary.json").read_text())["status"], "ANALYZED")

    def test_missing_structure_blocks_without_zero_scores(self):
        result = analyze_pair(self.sequence, self.ref, "A", self.root / "missing.pdb", self.ref, self.root / "blocked")
        self.assertEqual(result["status"], "BLOCKED")
        self.assertEqual(result["metrics"], {})

    def test_standard_metrics_require_explicit_reference_length(self):
        with patch("project1.metrics.run_tmscore") as tool:
            result = analyze_pair(self.sequence, self.ref, "A", self.ref, self.ref, self.root / "no_reference_length")
        tool.assert_not_called()
        self.assertEqual(result["status"], "ANALYZED")
        self.assertIsNone(result["reference_normalization_length"])
        self.assertIsNone(result["scored_fraction_of_reference_normalization"])
        self.assertIsNone(result["metrics"]["AUTO"]["tm_score"])
        self.assertIn("reference_length", result["metrics"]["AUTO"]["standard_tool"]["requirement"])

    def test_caught_scoring_error_cannot_leave_analyzed_success(self):
        with patch("project1.metrics.run_tmscore", side_effect=ValueError("SYNTHETIC scoring failure")):
            result = analyze_pair(self.sequence, self.ref, "A", self.ref, self.ref,
                                  self.root / "scoring_failure", reference_length=5)
        self.assertEqual(result["status"], "BLOCKED")
        self.assertEqual(result["metrics"], {})
        self.assertNotIn("paired_differences_DEEPMSA_minus_AUTO", result)
        self.assertIn("SYNTHETIC scoring failure", result["diagnostic"])


class TMscoreWrapperTests(unittest.TestCase):
    def test_documented_format_and_length_normalization(self):
        result = parse_tmscore_output(SYNTHETIC_TMSCORE_OUTPUT, reference_length=10, mask_count=5)
        self.assertEqual(result["tm_score"], 0.4567)
        self.assertEqual(result["gdt_ts"], 0.4)
        self.assertEqual(result["gdt_ts_native_mask_normalized"], 0.8)

    def test_wrong_normalization_or_lost_pairs_rejected(self):
        with self.assertRaises(ValueError):
            parse_tmscore_output(SYNTHETIC_TMSCORE_OUTPUT, 11, 5)
        with self.assertRaises(ValueError):
            parse_tmscore_output(SYNTHETIC_TMSCORE_OUTPUT.replace("in common=    5", "in common=    4"), 10, 5)
        with self.assertRaises(ValueError):
            parse_tmscore_output(SYNTHETIC_TMSCORE_OUTPUT.split("TM-score    = 0.4567")[0], 10, 5)

    def test_missing_executable_is_unavailable(self):
        with patch("project1.metrics.shutil.which", return_value=None):
            result = run_tmscore(Path("synthetic_model.pdb"), Path("synthetic_native.pdb"), 10, 5,
                                 "__nonexistent_official_TMscore_for_test__")
        self.assertEqual(result["status"], "UNAVAILABLE")
        self.assertIsNone(result["tm_score"])
        self.assertIsNone(result["gdt_ts"])
        self.assertIn("official", result["requirement"])

    @unittest.skipUnless(os.environ.get("PROJECT1_TMSCORE"), "Live official-tool equivalence unavailable: set PROJECT1_TMSCORE to an existing official TMscore executable; none installed by workflow")
    def test_live_wrapper_matches_direct_official_implementation(self):
        # Optional direct-reference validation; all coordinates remain synthetic.
        executable = os.environ["PROJECT1_TMSCORE"]
        with tempfile.TemporaryDirectory(prefix="project1_synthetic_TMscore_") as folder:
            root = Path(folder)
            coords = np.array([[i, np.sin(i), np.cos(i)] for i in range(20)])
            for name, xyz in (("native", coords), ("model", coords + [4, 6, 7])):
                _write_ca_pdb(root / f"{name}.pdb", list(range(1, 21)), "ACDEFGHIKLMNPQRSTVWYA", xyz)
            direct = subprocess.run([executable, str(root / "model.pdb"), str(root / "native.pdb"), "-l", "25"], capture_output=True, text=True, check=True)
            expected = parse_tmscore_output(direct.stdout, 25, 20)
            actual = run_tmscore(root / "model.pdb", root / "native.pdb", 25, 20, executable)
            self.assertEqual(actual["tm_score"], expected["tm_score"])
            self.assertEqual(actual["gdt_ts"], expected["gdt_ts"])


if __name__ == "__main__":
    unittest.main()
