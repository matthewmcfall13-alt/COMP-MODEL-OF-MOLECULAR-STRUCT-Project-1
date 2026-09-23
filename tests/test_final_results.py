"""Synthetic publishing fixtures only; never use the real scientific inputs."""
import csv
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from project1.final_results import CONDITIONS, publish_results


class FinalResultsTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="project1_SYNTHETIC_publication_")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.output = self.root / "outputs/analyses/uploaded/synthetic_hash"
        self.output.mkdir(parents=True)
        (self.root / "analysis").mkdir()
        self.table = self.root / "analysis/casp15_structure_comparison_metrics.csv"
        self.team_rows = []
        msas, structures = [], []
        for target, conditions in CONDITIONS.items():
            for index, condition in enumerate(conditions):
                self.team_rows.append({"Target": target, "Condition": condition, "PDB": "SYNTHETIC",
                                      "MSA_depth": 3 - min(index, 1), "pLDDT": 70 + index,
                                      "pTM": .7 + .01 * index, "TM_score": .8 + .02 * index,
                                      "lDDT_CA": 80 + index, "RMSD_CA_A": 2 - .25 * index,
                                      "RMSD100_A": 1 - .1 * index, "Mapped_residues": 5})
                msas.append({"target": target, "condition": condition.replace(" ", "_"),
                             "status": "ANALYZED", "sequence_count_including_query": 3 - min(index, 1),
                             "unique_query_coordinate_rows": 2, "distinct_non_query_sequences": 1,
                             "source_sha256": "a" * 64, "neff": None,
                             "neff_status": "UNAVAILABLE_EXACT_ROW_LIMIT"})
                structures.append({"target": target, "condition": condition.replace(" ", "_"),
                                   "target_length": 5, "scored_residue_count": 5, "scored_fraction_of_target": 1,
                                   "rmsd_ca_angstrom": 2 - .25 * index, "lddt_ca_common_mask": .8 + index / 100,
                                   "ca_bfactor_residue_mean_all_mapped": 70 + index})
        self.write_team_table()
        self.result = {"source_git_revision": "SYNTHETIC_revision", "generated_utc": "SYNTHETIC_time",
                       "status": "COMPLETED_FOR_AVAILABLE_UPLOADS", "originals_integrity": "UNCHANGED",
                       "issues": [], "alignments": msas, "structures": {"rows": structures},
                       "input_inventory": [{"path": str(self.table), "exists": True,
                                            "size_bytes": self.table.stat().st_size,
                                            "sha256": hashlib.sha256(self.table.read_bytes()).hexdigest()}]}
        (self.output / "a_plot.png").write_bytes(b"SYNTHETIC image bytes; copied, never rendered")
        (self.output / "query_coordinates.fasta").write_text(">synthetic\nACDEF\n")
        with (self.output / "derived_paths.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["source", "derived"])
            writer.writeheader()
            writer.writerow({"source": str(self.table), "derived": str(self.output / "a_plot.png")})
        (self.output / "input_inventory.csv").write_text("path\nPRIVATE_ABSOLUTE_PATH\n")

    def write_team_table(self):
        with self.table.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(self.team_rows[0]))
            writer.writeheader()
            writer.writerows(self.team_rows)

    @staticmethod
    def read_csv(path):
        with Path(path).open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))

    def test_complete_merge_signed_deltas_and_portable_derivatives_preserve_originals(self):
        before = self.table.read_bytes()
        dest = publish_results(self.root, self.result, self.output)
        self.assertEqual(dest, self.root / "results/matthew_analysis")
        combined = self.read_csv(dest / "combined_results.csv")
        self.assertEqual(len(combined), 7)
        self.assertTrue(all(row["availability"] == "COMPLETE" for row in combined))
        self.assertEqual(combined[1]["team_ptm"], "0.71")
        self.assertEqual(combined[1]["team_results_source"], "analysis/casp15_structure_comparison_metrics.csv")
        self.assertEqual(combined[1]["calculated_rmsd_ca_angstrom"], "1.75")
        self.assertEqual(combined[0]["neff"], "")
        pairs = self.read_csv(dest / "paired_custom_minus_auto.csv")
        self.assertEqual(len(pairs), 4)
        self.assertAlmostEqual(float(pairs[0]["team_tm_score_delta"]), .02)
        self.assertAlmostEqual(float(pairs[0]["team_rmsd_ca_angstrom_delta"]), -.25)
        self.assertEqual(pairs[1]["scope"], "Additional comparison")
        inventory = self.read_csv(dest / "source_inventory.csv")
        self.assertEqual(inventory[0]["path"], "analysis/casp15_structure_comparison_metrics.csv")
        self.assertEqual(inventory[0]["sha256"], hashlib.sha256(before).hexdigest())
        self.assertEqual(self.table.read_bytes(), before)
        self.assertEqual((dest / "a_plot.png").read_bytes(), (self.output / "a_plot.png").read_bytes())
        self.assertFalse((dest / "query_coordinates.fasta").exists())
        self.assertFalse((dest / "input_inventory.csv").exists())
        portable = self.read_csv(dest / "derived_paths.csv")[0]
        self.assertEqual(portable["source"], "../../analysis/casp15_structure_comparison_metrics.csv")
        self.assertEqual(portable["derived"], "./a_plot.png")
        report = (dest / "REPORT.md").read_text(encoding="utf-8")
        self.assertIn("complete for the seven supplied study conditions", report)
        self.assertNotIn(str(self.root), report)
        self.assertIn("accepted team results", report)

    def test_missing_condition_stays_partial_without_inventing_metrics(self):
        self.result["alignments"].pop()
        self.team_rows.pop()
        self.write_team_table()
        self.result["status"] = "PARTIAL_REQUIRES_REVIEW"
        self.result["issues"] = ["Synthetic missing condition"]
        dest = publish_results(self.root, self.result, self.output)
        missing = self.read_csv(dest / "combined_results.csv")[-1]
        self.assertEqual(missing["target"], "T1122")
        self.assertIn("ALIGNMENT", missing["availability"])
        self.assertIn("TEAM_RESULTS", missing["availability"])
        self.assertEqual(missing["team_tm_score"], "")
        self.assertEqual(missing["sequence_count_including_query"], "")
        pair = self.read_csv(dest / "paired_custom_minus_auto.csv")[-1]
        self.assertEqual(pair["availability"], "PARTIAL")
        self.assertEqual(pair["team_tm_score_delta"], "")
        report = (dest / "REPORT.md").read_text(encoding="utf-8")
        self.assertIn("**Status: partial", report)
        self.assertNotIn("complete for the seven", report)
        self.assertFalse(json.loads((dest / ".generated_manifest.json").read_text())["complete"])

    def test_refresh_removes_only_stale_owned_files_and_preserves_user_additions(self):
        dest = publish_results(self.root, self.result, self.output)
        user_note = dest / "allen_notes.txt"
        user_note.write_text("Keep these notes.")
        (self.output / "a_plot.png").unlink()
        publish_results(self.root, self.result, self.output)
        self.assertFalse((dest / "a_plot.png").exists())
        self.assertEqual(user_note.read_text(), "Keep these notes.")
        (dest / "REPORT.md").write_text("Manually edited: preserve me.")
        with self.assertRaisesRegex(ValueError, "was edited"):
            publish_results(self.root, self.result, self.output)
        self.assertEqual((dest / "REPORT.md").read_text(), "Manually edited: preserve me.")

    def test_nonempty_unowned_destination_is_not_overwritten(self):
        dest = self.root / "results/matthew_analysis"
        dest.mkdir(parents=True)
        colleague_file = dest / "REPORT.md"
        colleague_file.write_text("A colleague's existing result.")
        with self.assertRaisesRegex(ValueError, "ownership manifest"):
            publish_results(self.root, self.result, self.output)
        self.assertEqual(colleague_file.read_text(), "A colleague's existing result.")

    def test_duplicate_team_key_is_rejected_instead_of_arbitrary_selection(self):
        self.team_rows.append(dict(self.team_rows[0]))
        self.write_team_table()
        with self.assertRaisesRegex(ValueError, "Duplicate target/condition"):
            publish_results(self.root, self.result, self.output)


if __name__ == "__main__":
    unittest.main()
