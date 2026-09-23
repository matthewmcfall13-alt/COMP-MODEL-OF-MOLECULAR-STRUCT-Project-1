"""Current-evidence tables and figures. Never reads historical analysis folders."""
from pathlib import Path

from .inputs import sha256, write_json


def build_reports(root, state):
    from .workflow import write_csv
    root = Path(root)
    models, alignments, accuracy, availability = [], [], [], []
    analyzed, alignment_pairs, confidence_pairs = [], [], []
    for ident, row in state["targets"].items():
        alignment_by_condition, confidence_by_condition = {}, {}
        for condition in ("AUTO", "DEEPMSA"):
            stage = row[condition]
            source = row["deepmsa"] if condition == "DEEPMSA" else stage
            stats = source.get("alignment_analysis", {})
            if "query_length" in stats:
                alignment_by_condition[condition] = stats
                record = {"target": ident, "condition": condition,
                          **{k: stats.get(k) for k in ("query_length", "sequence_count_including_query", "gap_fraction", "mean_identity_to_query_excluding_query", "duplicate_fraction", "blosum62_normalized_sum_of_pairs", "blosum62_eligible_pair_denominator", "neff", "neff_status")}}
                alignments.append(record)
            for model in stage.get("import", {}).get("models", []):
                conf = model.get("confidence", {})
                primary = model["path"] == stage.get("import", {}).get("primary_model")
                models.append({"target": ident, "condition": condition, "model": model["path"], "primary": primary,
                               "verification_status": stage.get("import", {}).get("status"), "server_rank": model.get("server_rank"),
                               "ranking_score": model.get("ranking_score"), "ptm": conf.get("ptm"),
                               "atom_mean_plddt": conf.get("atom_mean_plddt"), "residue_mean_plddt": conf.get("residue_mean_plddt"),
                               "pae_mapping_status": conf.get("pae_mapping_status", "UNAVAILABLE")})
                if primary and stage.get("import", {}).get("comparison_eligible"):
                    confidence_by_condition[condition] = conf
                    fingerprint = sha256((model.get("sha256", "") + str(model.get("full_path", ""))).encode())[:16]
                    dest = root / "outputs" / "analyses" / ident / fingerprint / (condition + "_confidence")
                    dest.mkdir(parents=True, exist_ok=True)
                    write_json(dest / "confidence.json", conf)
                    state["review_artifacts"].append(str(dest))
                    profile = conf.get("residue_plddt", [])
                    if profile:
                        write_csv(dest / "residue_plddt.csv", profile, list(profile[0]))
                        from matplotlib.figure import Figure
                        fig = Figure(figsize=(8, 3.5), layout="constrained")
                        ax = fig.subplots()
                        ax.plot([int(r["label_seq_id"]) for r in profile], [r["plddt"] for r in profile])
                        ax.set(xlabel="Target residue position (mmCIF label_seq_id)", ylabel="Mean atom pLDDT within residue (0–100)", ylim=(0, 100), title=f"{ident} {condition} primary model confidence")
                        fig.savefig(dest / "plddt_profile.png", dpi=180)
                    if conf.get("pae") is not None and conf.get("pae_mapping_status") == "VERIFIED":
                        from matplotlib.figure import Figure
                        fig = Figure(figsize=(5.5, 4.5), layout="constrained")
                        ax = fig.subplots()
                        im = ax.imshow(conf["pae"], origin="lower", interpolation="nearest")
                        ax.set(xlabel="Token index (0-based; mapping in confidence.json)", ylabel="Token index (0-based)", title=f"{ident} {condition} PAE")
                        fig.colorbar(im, ax=ax, label="Predicted aligned error (Å)")
                        fig.savefig(dest / "pae.png", dpi=180)
        if len(alignment_by_condition) == 2:
            alignment_pairs.append(ident)
        measured = row.get("accuracy", {})
        if measured.get("status") == "ANALYZED":
            analyzed.append(ident)
            for condition, metric in measured.get("metrics", {}).items():
                accuracy.append({"target": ident, "condition": condition, "scored_residue_count": measured.get("scored_residue_count"),
                                 "scored_fraction_of_target": measured.get("scored_fraction_of_target"),
                                 **{k: metric.get(k) for k in ("rmsd_angstrom", "lddt_ca_common_mask", "tm_score", "gdt_ts", "molprobity")}})
            if len(confidence_by_condition) == 2 and any(all(c.get(key) is not None for c in confidence_by_condition.values()) for key in ("residue_mean_plddt", "ptm")):
                confidence_pairs.append(ident)
            plot_differences(ident, measured)
        # Each planned metric receives an explicit row even in the FASTA-only state.
        definitions = {
            "alignment_coverage": "Validated returned AUTO and DEEPMSA A3Ms",
            "identity_to_query": "Validated returned AUTO and DEEPMSA A3Ms",
            "sequence_redundancy": "Validated returned AUTO and DEEPMSA A3Ms",
            "normalized_blosum62_sp": "Validated returned AUTO and DEEPMSA A3Ms",
            "neff": "Validated A3Ms within the exact computation bound (500 rows)",
            "rmsd": "Validated paired primary models and confirmed experimental reference/mapping",
            "tm_score": "Official TMscore executable plus confirmed full reference length and paired inputs",
            "lddt_ca": "Validated paired primary models and experimental reference; common-mask C-alpha scope",
            "gdt_ts": "Official TMscore executable plus confirmed full reference length and paired inputs",
            "molprobity": "Separately supported MolProbity/Phenix geometry tool (not integrated)",
            "plddt": "Returned confidence with validated atom-to-residue mapping",
            "ptm": "Returned summary confidence pTM",
            "pae": "Returned PAE with verified token mapping",
        }
        metric_keys = {"rmsd": "rmsd_angstrom", "tm_score": "tm_score", "lddt_ca": "lddt_ca_common_mask", "gdt_ts": "gdt_ts", "molprobity": "molprobity"}
        for name, requirement in definitions.items():
            available = False
            if name in ("alignment_coverage", "identity_to_query", "sequence_redundancy", "normalized_blosum62_sp"):
                available = len(alignment_by_condition) == 2
                if name == "identity_to_query":
                    available = available and all(a.get("mean_identity_to_query_excluding_query") is not None for a in alignment_by_condition.values())
                if name == "normalized_blosum62_sp":
                    available = available and all(a.get("blosum62_normalized_sum_of_pairs") is not None and a.get("blosum62_eligible_pair_denominator", 0) > 0 for a in alignment_by_condition.values())
            elif name == "neff":
                available = len(alignment_by_condition) == 2 and all(a.get("neff_status") == "EXACT" for a in alignment_by_condition.values())
            elif name in metric_keys:
                available = measured.get("status") == "ANALYZED" and measured.get("metric_coverage", {}).get(metric_keys[name]) == "AVAILABLE"
            elif len(confidence_by_condition) == 2:
                if name == "plddt":
                    available = all(c.get("residue_mean_plddt") is not None for c in confidence_by_condition.values())
                elif name == "ptm":
                    available = all(c.get("ptm") is not None for c in confidence_by_condition.values())
                elif name == "pae":
                    available = all(c.get("pae_mapping_status") == "VERIFIED" for c in confidence_by_condition.values())
            display_names = {"normalized_blosum62_sp": "normalized_BLOSUM62_sum_of_pairs", "neff": "Neff", "rmsd": "CA_RMSD", "tm_score": "TM_score", "lddt_ca": "lDDT_CA", "gdt_ts": "GDT_TS", "molprobity": "MolProbity", "plddt": "residue_mean_pLDDT", "ptm": "pTM", "pae": "PAE"}
            availability.append({"target": ident, "metric": display_names.get(name, name), "status": "AVAILABLE" if available else "UNAVAILABLE",
                                 "evidence": "outputs/status.json" if available else "", "requirement": requirement,
                                 "note": "Availability requires both paired conditions; individual available outputs remain in their tables"})
    write_csv(root / "outputs" / "model_confidence.csv", models, ["target", "condition", "model", "primary", "verification_status", "server_rank", "ranking_score", "ptm", "atom_mean_plddt", "residue_mean_plddt", "pae_mapping_status"])
    write_csv(root / "outputs" / "alignment_summary.csv", alignments, ["target", "condition", "query_length", "sequence_count_including_query", "gap_fraction", "mean_identity_to_query_excluding_query", "duplicate_fraction", "blosum62_normalized_sum_of_pairs", "blosum62_eligible_pair_denominator", "neff", "neff_status"])
    write_csv(root / "outputs" / "accuracy_summary.csv", accuracy, ["target", "condition", "scored_residue_count", "scored_fraction_of_target", "rmsd_angstrom", "lddt_ca_common_mask", "tm_score", "gdt_ts", "molprobity"])
    state["metric_availability"] = availability
    write_csv(root / "outputs" / "metric_availability.csv", availability, ["target", "metric", "status", "evidence", "requirement", "note"])
    answers = {}
    for key, supported, evidence, answer in (
        ("accuracy", analyzed, "outputs/accuracy_summary.csv", "Paired accuracy metrics and mask coverage are reported per target and condition; inspect each metric's direction and missing-tool status."),
        ("alignments", alignment_pairs, "outputs/alignment_summary.csv", "Both returned alignments are available for these targets; compare depth, coverage, identity, redundancy and eligible-pair-normalized BLOSUM62."),
        ("confidence_accuracy", confidence_pairs, "outputs/model_confidence.csv", "Primary confidence and measured accuracy can be compared descriptively with accuracy_summary.csv; confidence is not accuracy and no inferential p-values are computed."),
    ):
        answers[key] = {"status": "PARTIALLY_ANSWERED" if supported else "UNANSWERED", "answer": (answer + " Available targets: " + ", ".join(supported)) if supported else "No complete real evidence for this question is available yet.",
                        "evidence": [evidence] if supported else [], "limitations": "Only available runs support conclusions. Missing references, unresolved controls, unavailable alignments and optional metrics remain explicit. Three targets do not establish general superiority."}
    answers["validation"] = {"status": "PARTIALLY_ANSWERED", "answer": "Input validation and request pairing are checked locally. Synthetic tests check software only; server acceptance and real paired mapping need downloaded evidence.",
                             "evidence": ["outputs/input_integrity.json"], "limitations": "Live optimized-score equivalence requires the official optional executable. Eligibility and experimental reference confirmation remain manual."}
    state["review_answers"] = answers


def plot_differences(ident, measured):
    from matplotlib.figure import Figure
    differences = measured.get("paired_differences_DEEPMSA_minus_AUTO", {})
    entries = [(key, value) for key, value in differences.items() if value is not None]
    if not entries:
        return
    fig = Figure(figsize=(max(5, 2.5 * len(entries)), 3.5), layout="constrained")
    axes = fig.subplots(1, len(entries), squeeze=False)[0]
    for ax, (name, value) in zip(axes, entries):
        ax.bar(["DEEPMSA − AUTO"], [value])
        ax.axhline(0, color="black", linewidth=0.8)
        ax.set_title(name)
        ax.set_ylabel("Difference (Å; lower better)" if name == "rmsd_angstrom" else "Difference (0–1 score; higher better)")
    fig.suptitle(f"{ident}: primary predictions, identical frozen residue mask")
    fig.savefig(Path(measured["analysis_directory"]) / "paired_accuracy_differences.png", dpi=180)
