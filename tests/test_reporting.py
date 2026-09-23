"""Synthetic reporting predicates: rejected evidence must not become study results."""
from project1.reporting import build_reports


def test_rejected_import_confidence_is_inventoried_not_compared(tmp_path):
    conf = {"ptm": 0.8, "residue_mean_plddt": 90, "pae": [[0]], "pae_mapping_status": "VERIFIED",
            "residue_plddt": [{"label_seq_id": "malformed", "plddt": 90}]}
    imp = {"status": "BLOCKED", "comparison_eligible": False, "primary_model": "synthetic.cif",
           "models": [{"path": "synthetic.cif", "confidence": conf}]}
    state = {"targets": {"SYNTHETIC": {"AUTO": {"import": imp}, "DEEPMSA": {"import": imp}, "deepmsa": {}}}, "review_artifacts": []}
    build_reports(tmp_path, state)
    assert all(r["status"] == "UNAVAILABLE" for r in state["metric_availability"] if r["metric"] in {"pTM", "PAE", "residue_mean_pLDDT"})
    assert len((tmp_path / "outputs/model_confidence.csv").read_text().splitlines()) == 3
    assert state["review_artifacts"] == []


def test_missing_alignment_scores_remain_unavailable(tmp_path):
    stats = {"query_length": 3, "mean_identity_to_query_excluding_query": None,
             "blosum62_normalized_sum_of_pairs": None, "blosum62_eligible_pair_denominator": 0}
    state = {"targets": {"SYNTHETIC": {"AUTO": {"alignment_analysis": stats}, "DEEPMSA": {},
                                       "deepmsa": {"alignment_analysis": stats}}}, "review_artifacts": []}
    build_reports(tmp_path, state)
    found = {r["metric"]: r["status"] for r in state["metric_availability"]}
    assert found["alignment_coverage"] == "AVAILABLE"
    assert found["identity_to_query"] == "UNAVAILABLE"
    assert found["normalized_BLOSUM62_sum_of_pairs"] == "UNAVAILABLE"
