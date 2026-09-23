"""Small rerunnable workflow: local preparation first, manual downloads later."""
from __future__ import annotations

import csv
import json
import traceback
from datetime import datetime, timezone
from pathlib import Path

from .inputs import (discover_fasta, file_hash, parse_a3m, select_a3m, server_request,
                     sha256, verify_pair, write_json, write_text)


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def load_config(config_path="config.json"):
    config_path = Path(config_path).resolve()
    config = json.loads(config_path.read_text(encoding="utf-8-sig"))
    root = Path(config.get("repo_root", ".")).expanduser()
    root = (config_path.parent / root).resolve() if not root.is_absolute() else root.resolve()
    data = Path(config.get("data_root") or "data").expanduser()
    data = (root / data).resolve() if not data.is_absolute() else data.resolve()
    config["_repo_root"], config["_data_root"] = root, data
    config["_config_path"] = config_path
    return config


def resolve_path(value, root):
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (root / path).resolve()


def write_csv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def prepare(config_path="config.json"):
    config = load_config(config_path)
    root, data = config["_repo_root"], config["_data_root"]
    records, issues = discover_fasta(root, config.get("input_overrides"), [data])
    selected = config.get("selected_targets", [])
    if not selected:
        issues.append("Select target IDs explicitly in config.json; no targets were chosen automatically")
    if len(selected) != len(set(selected)):
        raise ValueError("selected_targets contains duplicate IDs")
    if config.get("pilot") not in selected:
        issues.append("pilot must name a selected target")
    manifest_path = root / "target_manifest.json"
    old = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists() else {"targets": {}}
    manifest = old
    manifest.setdefault("targets", {})
    manifest["schema_version"] = 1
    manifest["candidate_records"] = [{k: v for k, v in r.items() if k != "sequence"} for r in records]
    manifest["discovery_issues"] = issues
    state = {"targets": {}, "issues": issues, "data_root": str(data), "repo_root": str(root),
             "evidence_type": "REAL_INPUTS_ONLY_NO_SERVER_RESULTS_ASSUMED", "generated_at": utc_now()}
    for ident in selected:
        candidates = [r for r in records if r["id"] == ident]
        previous = manifest["targets"].setdefault(ident, {})
        # Manual annotations live alongside regenerated integrity fields and survive reruns.
        for key, value in {"casp_eligibility": "UNRESOLVED", "reference_pdb": None, "reference_chain": None,
                           "reference_confirmed": False, "reference_assignment_evidence": None,
                           "manual": {}, "notes": ""}.items():
            previous.setdefault(key, value)
        previous["pilot"] = ident == config.get("pilot")
        previous["selected"] = True
        row = {"target_id": ident, "pilot": previous["pilot"], "AUTO": {"state": "BLOCKED"},
               "DEEPMSA": {"state": "WAITING_FOR_DEEPMSA"}, "deepmsa": {"state": "WAITING_FOR_DEEPMSA"}}
        state["targets"][ident] = row
        if len(candidates) != 1:
            row["errors"] = [f"Expected one FASTA record for {ident}; found {len(candidates)}"]
            previous["validation_status"] = "BLOCKED"
            continue
        rec = candidates[0]
        row["errors"] = rec["errors"]
        row["sequence"] = rec["sequence"]
        row["sequence_length"] = rec["sequence_length"]
        previous.update({k: v for k, v in rec.items() if k != "sequence"})
        previous["biology_flag"] = ("SUBUNIT_HEADER: confirm CASP single-chain eligibility; declared isolated-chain policy"
                                     if "subunit" in rec["header"].lower() else "Eligibility still requires CASP confirmation")
        row["biology_flag"] = previous["biology_flag"]
        if rec["errors"]:
            continue
        target_dir = root / "submissions" / ident
        fasta = target_dir / f"{ident}_DeepMSA2.fasta"
        sequence_file = target_dir / f"{ident}_sequence.txt"
        write_text(fasta, f">{ident}\n{rec['sequence']}\n")
        write_text(sequence_file, rec["sequence"] + "\n")
        request = server_request(ident, rec["sequence"], config.get("seed", 42))
        request_path = target_dir / (request[0]["name"] + ".json")
        write_json(request_path, request)
        row["AUTO"] = {"state": "PREPARED", "request_path": str(request_path),
                       "request_sha256": file_hash(request_path), "job_name": request[0]["name"]}
        row["deepmsa"].update({"fasta_path": str(fasta), "sequence_path": str(sequence_file),
                                  "proposed_mode": config.get("proposed_deepmsa_mode", "fast")})
        settings = config.get("targets", {}).get(ident, {})
        for stage in ("deepmsa", "AUTO", "DEEPMSA"):
            # Windows is case-insensitive: deepmsa and DEEPMSA cannot name distinct folders.
            stage_path = Path("deepmsa2") if stage == "deepmsa" else Path("alphafold") / stage
            default_dir = data / "downloads" / ident / stage_path / "run1"
            default_dir.mkdir(parents=True, exist_ok=True)
            supplied = settings.get(stage, {})
            download = resolve_path(supplied.get("input_path", str(default_dir)), root)
            row[stage]["download_dir"] = str(download)
            # Actual mode, IDs, timestamps are user-entered evidence, not inferred from a generated file.
            manual = previous["manual"].setdefault(stage, {})
            for key in ("submitted_at", "job_id", "actual_mode", "returned_seed", "notes", "submitted_request_sha256", "submitted_sequence_sha256"):
                if key in supplied:
                    manual[key] = supplied[key]
            row[stage]["manual"] = dict(manual)
            if manual.get("submitted_at"):
                submitted_hash = supplied.get("submitted_request_sha256", manual.get("submitted_request_sha256"))
                if stage == "deepmsa":
                    submitted_hash = supplied.get("submitted_sequence_sha256", manual.get("submitted_sequence_sha256"))
                    current_hash = rec["sequence_sha256"]
                    hash_key = "submitted_sequence_sha256"
                else:
                    current_hash = row[stage].get("request_sha256")
                    hash_key = "submitted_request_sha256"
                if submitted_hash is None and current_hash:
                    manual[hash_key] = current_hash
                    submitted_hash = current_hash
                if current_hash and submitted_hash == current_hash:
                    row[stage]["state"] = "SUBMITTED"
                elif current_hash:
                    row[stage]["submission_warning"] = "Previous submission metadata belongs to a different input hash"
        previous["prepared_auto_sha256"] = row["AUTO"]["request_sha256"]
    # Deselection is recorded without discarding annotations.
    for ident, entry in manifest["targets"].items():
        entry["selected"] = ident in selected
    write_json(manifest_path, manifest)
    write_json(root / "outputs" / "status.json", state)
    write_json(root / "outputs" / "input_integrity.json", {
        "files": [{"path": r["source_path"], "sha256": r["source_sha256"],
                   "sequence_sha256": r["sequence_sha256"], "errors": r["errors"]} for r in records],
        "issues": issues, "originals_modified": False})
    write_start_here(root, config, state)
    print_status(state, config.get("pilot"))
    return config, state


def write_start_here(root, config, state):
    pilot = config.get("pilot")
    row = state["targets"].get(pilot, {})
    d = row.get("deepmsa", {})
    action, action_path, action_destination = next_action(state, pilot)
    lines = ["# Start here — manual server workflow", "", f"Pilot: **{pilot}**. {config.get('pilot_rationale', '')}", "",
             "## Current next action", "", action, f"Input/path: `{action_path}`", f"Destination: `{action_destination}`", "",
             "## First manual action", "", "Open https://zhanggroup.org/DeepMSA2/ yourself.",
             f"Upload `{d.get('fasta_path', 'BLOCKED: resolve input selection')}` if the UI offers a file input.",
             f"If it accepts pasted sequence only, paste the complete contents of `{d.get('sequence_path', 'BLOCKED')}`.",
             f"Proposed mode: **{config.get('proposed_deepmsa_mode', 'fast')}** for every target. The official page returned HTTP 403 during preparation on 2026-09-22; verify the current mode names, accepted length, quotas and submission controls in the UI. Do not truncate or change the sequence if rejected.",
             "Enter only this target, complete any required fields yourself, and record actual mode, job ID and submission timestamp in config.json (README has examples). Do not put passwords or access-bearing result URLs in config.",
             f"Download the complete real result archive/alignment to `{d.get('download_dir', 'BLOCKED')}`. Keep its original bytes and filename.",
             "Rerun **Project1_Server_Workflow.ipynb → Restart Kernel and Run All**. Multiple A3Ms require an explicit selection or a documented service-final designation in config.", "",
             "## Pilot AlphaFold Server pair", "", "Open https://alphafoldserver.com/ yourself and use its JSON upload control.",
             f"AUTO upload: `{row.get('AUTO', {}).get('request_path', 'BLOCKED')}`.",
             f"DEEPMSA upload: `{row.get('DEEPMSA', {}).get('request_path', 'WAITING_FOR_DEEPMSA — no placeholder JSON is created')}`.",
             f"Both requests use one protein copy, no additional entities, templates disabled, requested seed {config.get('seed', 42)}. AUTO omits unpairedMsa; DEEPMSA embeds the validated full A3M. If the UI rejects a setting, record the rejection; do not silently alter the comparison.",
             f"Save the **complete AUTO ZIP** to `{row.get('AUTO', {}).get('download_dir', 'BLOCKED')}`.",
             f"Save the **complete DEEPMSA ZIP** to `{row.get('DEEPMSA', {}).get('download_dir', 'BLOCKED')}`.",
             "Record submissions separately from preparation. Local validation has not established website acceptance or a successful job.", "",
             "## All selected targets", "", "See `outputs/manual_submission_sheet.csv` for exact request names, paths, settings and current states.",
             "Reference mappings and CASP eligibility remain UNRESOLVED. Headers identify T1106s1 and T1151s2 as complex subunits; submissions use the declared isolated single-chain policy. Confirm their assignment eligibility before interpreting them as eligible single-chain CASP targets.", "",
             "## Run locally", "", "From this repository: `.\\.venv\\Scripts\\python.exe -m project1 --config config.json`.",
             "The notebook and command perform local validation/import/analysis only. No websites are submitted or polled. Missing downloads are ordinary waiting states.", ""]
    write_text(root / "START_HERE.md", "\n".join(lines))
    sheet = []
    for ident, target in state["targets"].items():
        for condition in ("deepmsa", "AUTO", "DEEPMSA"):
            stage = target[condition]
            sheet.append({"target": ident, "pilot": target["pilot"], "condition": condition,
                          "state": stage["state"], "input_file": stage.get("request_path", stage.get("fasta_path", "")),
                          "download_directory": stage.get("download_dir", ""), "job_name": stage.get("job_name", ""),
                          "seed": config.get("seed", 42) if condition != "deepmsa" else "",
                          "templates": "disabled" if condition != "deepmsa" else "n/a",
                          "chain_copies": 1, "mode": stage.get("manual", {}).get("actual_mode", config.get('proposed_deepmsa_mode', 'fast') + " proposed; UI verification required") if condition == "deepmsa" else "",
                          "errors": "; ".join(stage.get("issues", []) + target.get("errors", []))})
    write_csv(root / "outputs" / "manual_submission_sheet.csv", sheet,
              ["target", "pilot", "condition", "state", "input_file", "download_directory", "job_name", "seed", "templates", "chain_copies", "mode", "errors"])


def next_action(state, pilot):
    p = state["targets"].get(pilot, {})
    deep = p.get("deepmsa", {})
    if not p or p.get("errors"):
        return "Resolve blocked FASTA inputs in config.json", "config.json", "local notebook"
    if deep.get("state") == "BLOCKED":
        return "Resolve the DeepMSA import/selection issue in outputs/status.json", deep.get("download_dir", ""), "config.json"
    if deep.get("state") != "VALIDATED":
        if deep.get("state") in ("SUBMITTED", "IMPORTED"):
            return "Save the complete DeepMSA result and rerun the notebook", deep.get("download_dir", ""), "local notebook"
        return "Submit this single target to DeepMSA2 manually", deep.get("fasta_path", ""), "https://zhanggroup.org/DeepMSA2/"
    for condition in ("AUTO", "DEEPMSA"):
        stage = p.get(condition, {})
        if stage.get("state") == "PREPARED":
            return f"Upload the pilot {condition} JSON manually", stage.get("request_path", ""), "https://alphafoldserver.com/"
        if stage.get("state") == "SUBMITTED":
            return f"Save the complete pilot {condition} result ZIP and rerun", stage.get("download_dir", ""), "local notebook"
        if stage.get("state") in ("BLOCKED", "IMPORTED_NEEDS_REVIEW", "IMPORTED_UNVERIFIED"):
            return f"Resolve pilot {condition} verification/ranking in outputs/status.json", stage.get("download_dir", ""), "config.json"
    return "Resolve reference/eligibility metadata, review local results, then continue the other prepared targets", "outputs/reference_resolution.csv", "config.json"


def print_status(state, pilot):
    print("\nTarget       AUTO                 DEEPMSA              DeepMSA2")
    for ident, row in state["targets"].items():
        print(f"{ident:12} {row['AUTO']['state']:20} {row['DEEPMSA']['state']:20} {row['deepmsa']['state']}")
        for issue in row.get("errors", []):
            print("  BLOCKED:", issue)
    p = state["targets"].get(pilot, {})
    action, path, destination = next_action(state, pilot)
    print("\nNext action:", action)
    print("Input/path:", path)
    print("Destination:", destination)
    print("Save DeepMSA2 downloads in:", p.get("deepmsa", {}).get("download_dir", "unresolved"))
    print("AUTO JSON:", p.get("AUTO", {}).get("request_path", "blocked"))
    if p.get("DEEPMSA", {}).get("request_path"):
        print("Validated DEEPMSA JSON:", p["DEEPMSA"]["request_path"])
    if destination == "https://alphafoldserver.com/":
        stage = next((p[c] for c in ("AUTO", "DEEPMSA") if p.get(c, {}).get("request_path") == path), {})
        print("Save this AlphaFold result ZIP in:", stage.get("download_dir", "unresolved"))


def run(config_path="config.json", package=True):
    config, state = prepare(config_path)
    root, data = config["_repo_root"], config["_data_root"]
    from .importers import ingest_files, import_alphafold
    from .metrics import alignment_statistics, analyze_pair
    # Every result analysis gets an input-derived directory; current.json points only to this run.
    current_artifacts, reference_rows = [], []
    manifest_path = root / "target_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for ident, row in state["targets"].items():
        if row.get("errors"):
            continue
        options = config.get("targets", {}).get(ident, {})
        try:
            ds = options.get("deepmsa", {})
            incoming = Path(row["deepmsa"]["download_dir"])
            imported = ingest_files(incoming, data / "retained" / ident / "deepmsa2")
            row["deepmsa"]["import"] = imported
            candidate_entries = [f for f in imported.get("files", []) if Path(f["path"]).suffix.lower() == ".a3m"]
            # An archive plus its already-extracted identical file is one alignment candidate.
            candidates = [Path(f["path"]) for f in {f["sha256"]: f for f in candidate_entries}.values()]
            row["deepmsa"]["candidate_inventory"] = [str(p) for p in candidates]
            if imported.get("status") == "BLOCKED":
                raise ValueError("DeepMSA archive import blocked: " + "; ".join(imported.get("issues", [])))
            chosen, rationale = select_a3m(candidates, ds.get("a3m_selection"), ds.get("service_final"))
            if chosen:
                parsed = parse_a3m(chosen.read_bytes().decode("utf-8"), row["sequence"])
                row["deepmsa"].update({"state": "VALIDATED", "chosen_alignment": str(chosen),
                                       "alignment_sha256": file_hash(chosen), "selection_rationale": rationale,
                                       "search_metadata": ds.get("search_metadata", {}),
                                       "actual_mode": row["deepmsa"]["manual"].get("actual_mode", "UNVERIFIED")})
                request = server_request(ident, row["sequence"], config.get("seed", 42), parsed["text"])
                auto = json.loads(Path(row["AUTO"]["request_path"]).read_text(encoding="utf-8"))
                verify_pair(auto, request)
                request_path = root / "submissions" / ident / (request[0]["name"] + ".json")
                write_json(request_path, request)
                row["DEEPMSA"].update({"state": "PREPARED", "request_path": str(request_path),
                                        "request_sha256": file_hash(request_path), "job_name": request[0]["name"], "paired_request_check": "VALIDATED"})
                manual = manifest["targets"][ident]["manual"].setdefault("DEEPMSA", {})
                if manual.get("submitted_at"):
                    submitted_hash = options.get("DEEPMSA", {}).get("submitted_request_sha256", manual.get("submitted_request_sha256"))
                    if submitted_hash is None:
                        submitted_hash = row["DEEPMSA"]["request_sha256"]
                        manual["submitted_request_sha256"] = submitted_hash
                    if submitted_hash == row["DEEPMSA"]["request_sha256"]:
                        row["DEEPMSA"]["state"] = "SUBMITTED"
                    else:
                        row["DEEPMSA"]["submission_warning"] = "Prior submitted request differs from current validated A3M request"
                row["DEEPMSA"]["manual"] = dict(manual)
                analysis_dir = root / "outputs" / "analyses" / ident / parsed["sha256"][:16] / "DEEPMSA_alignment"
                write_text(analysis_dir / "query_coordinates.fasta", "".join(f">{h}\n{s}\n" for h, s in zip(parsed["headers"], parsed["analysis_rows"])))
                row["deepmsa"]["alignment_analysis"] = alignment_statistics(parsed["analysis_rows"], output_dir=analysis_dir)
                row["deepmsa"]["analysis_directory"] = str(analysis_dir)
                current_artifacts.extend([str(chosen), str(analysis_dir)])
            elif imported.get("files"):
                row["deepmsa"].update({"state": "IMPORTED", "issues": ["No A3M found in the returned files; select or obtain the real A3M output"]})
        except (ValueError, OSError, UnicodeError) as exc:
            row["deepmsa"].update({"state": "BLOCKED", "issues": [str(exc)]})
            row["DEEPMSA"]["state"] = "BLOCKED"
        except Exception as exc:
            diagnostic(root, ident + "_deepmsa", exc)
            row["deepmsa"].update({"state": "BLOCKED", "issues": [f"Unexpected {type(exc).__name__}: {exc}; see outputs/diagnostics"]})
        for condition in ("AUTO", "DEEPMSA"):
            stage = row[condition]
            if not stage.get("request_path"):
                continue
            try:
                expected = json.loads(Path(stage["request_path"]).read_text(encoding="utf-8"))
                imported = import_alphafold(Path(stage["download_dir"]), data / "retained" / ident / "alphafold" / condition,
                                             expected, options.get(condition, {}).get("primary_override"))
                stage["import"] = imported
                if imported.get("models") or imported.get("originals") or imported.get("files"):
                    stage["state"] = imported.get("status", "IMPORTED")
                if imported.get("primary_model"):
                    current_artifacts.append(imported["primary_model"])
                    primary = next(m for m in imported["models"] if m["path"] == imported["primary_model"])
                    current_artifacts.extend(primary[k] for k in ("summary_path", "full_path") if primary.get(k))
                current_artifacts.extend(r["path"] for r in imported.get("returned_requests", []))
                if condition == "AUTO":
                    try:
                        automatic_alignment(row, imported, root, current_artifacts, options.get("AUTO", {}))
                    except (ValueError, OSError, UnicodeError) as exc:
                        stage["alignment_analysis"] = {"status": "UNAVAILABLE", "reason": str(exc)}
            except (ValueError, OSError) as exc:
                stage.update({"state": "BLOCKED", "issues": [str(exc)]})
            except Exception as exc:
                diagnostic(root, ident + "_" + condition, exc)
                stage.update({"state": "BLOCKED", "issues": [f"Unexpected {type(exc).__name__}: {exc}; see outputs/diagnostics"]})
        prior = manifest["targets"][ident]
        ref = dict(prior.get("reference", {}))
        for key, field in (("pdb_id", "reference_pdb"), ("chain", "reference_chain"), ("confirmed", "reference_confirmed"),
                           ("assignment_evidence", "reference_assignment_evidence"), ("experimental", "reference_experimental"),
                           ("path", "reference_path"), ("reference_length", "reference_length")):
            if field in prior and prior[field] is not None:
                ref.setdefault(key, prior[field])
        ref.update(options.get("reference", {}))
        reference_rows.append({"target": ident, "pdb_id": ref.get("pdb_id", "UNRESOLVED"),
                               "chain": ref.get("chain", "UNRESOLVED"), "confirmed": ref.get("confirmed", False),
                               "assignment_evidence": ref.get("assignment_evidence", "UNRESOLVED"), "path": ref.get("path", "")})
        row["accuracy"] = {"status": "WAITING_FOR_CONFIRMED_REFERENCE_AND_VALIDATED_PAIR"}
        imports = [row[c].get("import", {}) for c in ("AUTO", "DEEPMSA")]
        if ref.get("confirmed") and ref.get("experimental") and ref.get("assignment_evidence") and ref.get("path") and ref.get("chain"):
            if all(i.get("primary_model") and i.get("comparison_eligible") for i in imports):
                try:
                    reference = resolve_path(ref["path"], root)
                    fingerprint = sha256("|".join([file_hash(reference)] + [file_hash(Path(i["primary_model"])) for i in imports] +
                                         [json.dumps(ref, sort_keys=True), json.dumps(options, sort_keys=True), str(config.get("tmscore_executable"))]).encode())[:16]
                    output_dir = root / "outputs" / "analyses" / ident / fingerprint / "accuracy"
                    row["accuracy"] = analyze_pair(row["sequence"], reference, ref["chain"],
                        imports[0]["primary_model"], imports[1]["primary_model"], output_dir,
                        auto_chain=options.get("AUTO", {}).get("chain", "A"), deep_chain=options.get("DEEPMSA", {}).get("chain", "A"),
                        reference_mapping=ref.get("mapping"), auto_mapping=options.get("AUTO", {}).get("mapping"),
                        deep_mapping=options.get("DEEPMSA", {}).get("mapping"), reference_length=ref.get("reference_length"),
                        tmscore_executable=config.get("tmscore_executable"), reference_model_index=ref.get("model_index"))
                    row["accuracy"]["analysis_directory"] = str(output_dir)
                    current_artifacts.extend([str(reference), str(output_dir)])
                except Exception as exc:
                    diagnostic(root, ident + "_accuracy", exc)
                    row["accuracy"] = {"status": "BLOCKED", "issues": [f"{type(exc).__name__}: {exc}"]}
    write_csv(root / "outputs" / "reference_resolution.csv", reference_rows,
              ["target", "pdb_id", "chain", "confirmed", "assignment_evidence", "path"])
    state["review_artifacts"] = current_artifacts
    write_json(manifest_path, manifest)
    from .reporting import build_reports
    build_reports(root, state)
    write_json(root / "outputs" / "current.json", {"generated_at": utc_now(), "artifacts": current_artifacts,
                                                      "note": "Only these artifacts belong to current validated inputs; older analyses remain historical"})
    write_json(root / "outputs" / "status.json", state)
    write_start_here(root, config, state)
    print_status(state, config.get("pilot"))
    if package:
        from .review import create_review
        state["review_zip"] = str(create_review(root, data, state))
        print("Review ZIP:", state["review_zip"])
    return state


def automatic_alignment(row, imported, root, current_artifacts, settings):
    candidates = imported.get("returned_alignments", [])
    candidates = [Path(x["path"] if isinstance(x, dict) else x) for x in candidates
                  if Path(x["path"] if isinstance(x, dict) else x).suffix.lower() == ".a3m"]
    if not candidates:
        row["AUTO"]["alignment_analysis"] = {"status": "UNAVAILABLE", "reason": "Automatic alignment was not returned"}
        return
    chosen, rationale = select_a3m(candidates, settings.get("a3m_selection"))
    parsed = parse_a3m(chosen.read_bytes().decode("utf-8"), row["sequence"])
    from .metrics import alignment_statistics
    output_dir = root / "outputs" / "analyses" / row["target_id"] / parsed["sha256"][:16] / "AUTO_alignment"
    row["AUTO"]["alignment_analysis"] = alignment_statistics(parsed["analysis_rows"], output_dir=output_dir)
    row["AUTO"]["alignment_selection_rationale"] = rationale
    row["AUTO"]["analysis_directory"] = str(output_dir)
    current_artifacts.extend([str(chosen), str(output_dir)])


def diagnostic(root, stage, exception):
    write_text(root / "outputs" / "diagnostics" / f"{stage}.txt", traceback.format_exc())
    print(f"BLOCKED {stage}: {type(exception).__name__}: {exception}")
