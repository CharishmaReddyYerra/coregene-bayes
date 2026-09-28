"""Command-line interface."""

from __future__ import annotations

import argparse
from pathlib import Path

from coregene_bayes.config import load_config, project_root, resolve_path
from coregene_bayes.data_sources import build_real_snapshot
from coregene_bayes.pipeline import run_pipeline


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="coregene-bayes")
    subparsers = parser.add_subparsers(dest="command", required=True)
    run = subparsers.add_parser("run", help="run the frozen real-data analysis")
    run.add_argument("--config", type=Path)
    run.add_argument("--data", type=Path)
    run.add_argument("--output", type=Path)

    fetch = subparsers.add_parser("fetch", help="download and compact all public inputs")
    fetch.add_argument("--config", type=Path)
    fetch.add_argument("--output", type=Path)
    fetch.add_argument("--cache", type=Path, default=Path(".cache/coregene-bayes"))
    fetch.add_argument("--insecure", action="store_true")
    fetch.add_argument("--existing-liver-snapshot", type=Path)
    return parser


def main() -> None:
    args = _parser().parse_args()
    if args.command == "run":
        outputs = run_pipeline(
            config_path=args.config,
            data_directory=args.data,
            output_directory=args.output,
        )
        print(f"Completed analysis; report: {outputs['report.html']}")
        return

    root = project_root()
    config = load_config(args.config)
    output = args.output or resolve_path(root, config["data_directory"])
    cache = args.cache if args.cache.is_absolute() else (root / args.cache).resolve()
    bayesv2g = resolve_path(root, config["bayesv2g_directory"])
    existing = args.existing_liver_snapshot
    if existing is None:
        candidate = root.parent / "bayesv2g" / "data" / "snapshots" / "sort1_case" / "eqtl_catalogue_qtd000266_liver_sort1_genes.tsv.gz"
        existing = candidate if candidate.exists() else None
    result = build_real_snapshot(
        output=output,
        cache=cache,
        endpoint=config["open_targets_endpoint"],
        l2g_threshold=float(config["l2g_candidate_threshold"]),
        p_threshold=float(config["candidate_p_threshold"]),
        bayesv2g_directory=bayesv2g,
        existing_liver_snapshot=existing,
        verify_tls=not args.insecure,
    )
    print(f"Wrote {result['candidate_selection']['unique_candidate_genes']} candidates to {output}")

