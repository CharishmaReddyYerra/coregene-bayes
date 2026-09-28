#!/usr/bin/env python3
"""Fetch and compact the complete public inputs for the real-data analysis."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from coregene_bayes.config import load_config, resolve_path  # noqa: E402
from coregene_bayes.data_sources import build_real_snapshot  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--cache", type=Path, default=Path(".cache/coregene-bayes"))
    parser.add_argument("--insecure", action="store_true")
    parser.add_argument("--existing-liver-snapshot", type=Path)
    args = parser.parse_args()
    config = load_config(args.config)
    output = args.output or resolve_path(ROOT, config["data_directory"])
    cache = args.cache if args.cache.is_absolute() else (ROOT / args.cache).resolve()
    bayesv2g = resolve_path(ROOT, config["bayesv2g_directory"])
    existing = args.existing_liver_snapshot
    if existing is None:
        candidate = ROOT.parent / "bayesv2g" / "data" / "snapshots" / "sort1_case" / "eqtl_catalogue_qtd000266_liver_sort1_genes.tsv.gz"
        existing = candidate if candidate.exists() else None
    manifest = build_real_snapshot(
        output=output,
        cache=cache,
        endpoint=config["open_targets_endpoint"],
        l2g_threshold=float(config["l2g_candidate_threshold"]),
        p_threshold=float(config["candidate_p_threshold"]),
        bayesv2g_directory=bayesv2g,
        existing_liver_snapshot=existing,
        verify_tls=not args.insecure,
    )
    print(f"Wrote real-data snapshot to {output}")
    print(manifest["counts"])


if __name__ == "__main__":
    main()

