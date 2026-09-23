"""Synthetic A3M fixtures in temporary directories, never study results."""
import csv
import hashlib
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from project1.uploaded_msa import (
    analyze_uploaded_msa, canonical_conservation_profile, read_uploaded_a3m,
)


class UploadedMSATests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="project1_synthetic_uploaded_msa_")
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def write(self, text):
        path = self.root / "input_SYNTHETIC.a3m"
        path.write_text(text, encoding="utf-8")
        return path

    def analyze(self, text, query="ACD"):
        source = self.write(text)
        with patch("project1.uploaded_msa._figures", return_value={}):
            return analyze_uploaded_msa(source, query, self.root / "derived")

    def test_colabfold_header_insertions_and_dots(self):
        source = self.write("#3\t1\n>query\nACD\n>h1 description\nAr.C-D\n")
        # Removing one insertion and one dot leaves four columns: rejected.
        with self.assertRaisesRegex(ValueError, "width"):
            read_uploaded_a3m(source, "ACD")
        source = self.write("#3\t1\n>query\nACD\n>h1 description\nAr.C-\n")
        names, rows, meta = read_uploaded_a3m(source, "ACD")
        self.assertEqual(names, ["query", "h1 description"])
        self.assertEqual(rows, ["ACD", "AC-"])
        self.assertEqual(meta["colabfold_metadata"], {"chain_lengths": [3], "copy_counts": [1]})
        self.assertEqual(meta["lowercase_insertion_residues_removed"], 1)
        self.assertEqual(meta["dot_placeholders_removed"], 1)

    def test_reject_incorrect_and_multimer_headers(self):
        for header in ("#4\t1", "#3\t2", "#3,3\t1,1", "#3\t1,1", "#something"):
            with self.subTest(header=header):
                source = self.write(header + "\n>query\nACD\n")
                with self.assertRaisesRegex(ValueError, "single-chain"):
                    read_uploaded_a3m(source, "ACD")

    def test_reject_query_mismatch_and_invalid_records(self):
        for text in (">query\nACE\n", ">query\nACD\n>empty\n", "ACD\n", ">query\nACD\n>bad\nAC*\n", ">query\nACD\n#3\t1\n", ">\nACD\n"):
            with self.subTest(text=text):
                with self.assertRaises(ValueError):
                    read_uploaded_a3m(self.write(text), "ACD")

    def test_query_only_is_real_analysis_but_no_homolog_evidence(self):
        summary = self.analyze(">query\nACD\n")
        self.assertEqual(summary["status"], "QUERY_ONLY")
        self.assertEqual(summary["homolog_row_count"], 0)
        self.assertEqual(summary["distinct_non_query_sequences"], 0)
        self.assertIsNone(summary["blosum62_normalized_sum_of_pairs"])
        self.assertEqual(summary["neff"], 1)
        with (self.root / "derived/alignment_positions.csv").open() as handle:
            positions = list(csv.DictReader(handle))
        self.assertEqual(positions[0]["non_query_consensus_fraction"], "")
        self.assertEqual(positions[0]["non_query_canonical_count"], "0")

    def test_duplicate_queries_do_not_imply_homolog_diversity(self):
        summary = self.analyze(">query\nACD\n>dup1\nACD\n>dup2\nACD\n")
        self.assertEqual(summary["status"], "NO_NON_QUERY_DIVERSITY")
        self.assertEqual(summary["sequence_count_including_query"], 3)
        self.assertEqual(summary["distinct_non_query_sequences"], 0)
        self.assertEqual(summary["duplicate_query_rows"], 2)
        self.assertEqual(summary["unique_query_coordinate_rows"], 1)
        self.assertEqual(summary["neff"], 1)

    def test_distinct_count_and_full_input_statistics(self):
        summary = self.analyze(">query\nACD\n>h1\nADD\n>h2\nADD\n>dup\nACD\n")
        self.assertEqual(summary["status"], "ANALYZED")
        self.assertEqual(summary["distinct_non_query_sequences"], 1)
        self.assertEqual(summary["duplicate_query_rows"], 1)
        self.assertEqual(summary["duplicate_rows_beyond_first"], 2)
        self.assertAlmostEqual(summary["mean_identity_to_query_excluding_query"], 7 / 9)
        self.assertNotIn("positions", summary)
        self.assertNotIn("rows", summary)

    def test_entropy_denominator_excludes_gaps_and_ambiguous(self):
        profile = canonical_conservation_profile(["A--", "C-X", "---", "X--"])
        self.assertEqual(profile[0]["canonical_count"], 2)
        self.assertAlmostEqual(profile[0]["shannon_entropy_bits"], 1)
        self.assertAlmostEqual(profile[0]["shannon_entropy_normalized_log2_20"], 1 / math.log2(20))
        self.assertEqual(profile[0]["consensus_fraction"], 0.5)
        for index in (1, 2):
            self.assertEqual(profile[index]["canonical_count"], 0)
            self.assertIsNone(profile[index]["shannon_entropy_bits"])
            self.assertIsNone(profile[index]["conservation_one_minus_normalized_entropy"])
            self.assertIsNone(profile[index]["consensus_fraction"])

    def test_uniform_20_amino_acids_have_maximum_entropy(self):
        profile = canonical_conservation_profile(list("ACDEFGHIKLMNPQRSTVWY"))[0]
        self.assertAlmostEqual(profile["shannon_entropy_normalized_log2_20"], 1)
        self.assertAlmostEqual(profile["conservation_one_minus_normalized_entropy"], 0)
        self.assertAlmostEqual(profile["consensus_fraction"], 1 / 20)

    def test_exact_neff_is_missing_over_500_rows(self):
        summary = self.analyze("".join(f">synthetic_{i}\nACD\n" for i in range(501)))
        self.assertEqual(summary["sequence_count_including_query"], 501)
        self.assertIsNone(summary["neff"])
        self.assertEqual(summary["neff_status"], "UNAVAILABLE_EXACT_ROW_LIMIT")
        self.assertFalse(summary["neff_approximation"])
        self.assertEqual(summary["duplicate_query_rows"], 500)

    def test_original_preserved_and_real_figures_and_artifacts_written(self):
        source = self.write("#3\t1\n>query\nACD\n>h1\nAq.-D\n")
        original = source.read_bytes()
        summary = analyze_uploaded_msa(source, "ACD", self.root / "derived")
        self.assertEqual(source.read_bytes(), original)
        self.assertEqual(summary["source_sha256"], hashlib.sha256(original).hexdigest())
        self.assertEqual((self.root / "derived/query_coordinates.fasta").read_text(), ">query\nACD\n>h1\nA-D\n")
        self.assertEqual(summary["plot_status"], "AVAILABLE")
        self.assertEqual(summary["overview_row_indices_zero_based"], [0, 1])
        for name, path in summary["artifacts"].items():
            self.assertTrue(Path(path).is_file(), name)
            self.assertGreater(Path(path).stat().st_size, 0, name)
        saved = json.loads((self.root / "derived/alignment_summary.json").read_text())
        self.assertEqual(saved, summary)

    def test_input_cannot_be_overwritten_by_output(self):
        source = self.root / "query_coordinates.fasta"
        source.write_text(">query\nACD\n")
        with self.assertRaisesRegex(ValueError, "overwrite"):
            analyze_uploaded_msa(source, "ACD", self.root)


if __name__ == "__main__":
    unittest.main()
