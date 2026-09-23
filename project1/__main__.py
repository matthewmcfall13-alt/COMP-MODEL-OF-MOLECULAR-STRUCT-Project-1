import argparse
from .workflow import prepare, run

parser = argparse.ArgumentParser(description="Local AF2 project analysis and retained manual-server preparation; never submits jobs")
parser.add_argument("--config", default="config.json")
parser.add_argument("--prepare-only", action="store_true")
parser.add_argument("--no-package", action="store_true")
parser.add_argument("--mode", choices=("auto", "uploaded", "manual-server"), default="auto",
                    help="auto audits configured uploaded results when present; otherwise prepares manual-server inputs")
args = parser.parse_args()
if args.prepare_only:
    prepare(args.config)
else:
    from .uploaded_analysis import uploads_available, run_uploaded, ensure_uploaded_source
    if args.mode == "auto":
        ensure_uploaded_source(args.config)
    if args.mode == "uploaded" or (args.mode == "auto" and uploads_available(args.config)):
        run_uploaded(args.config, package=not args.no_package)
    else:
        run(args.config, package=not args.no_package)
