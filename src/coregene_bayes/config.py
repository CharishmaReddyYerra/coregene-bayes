"""Configuration loading and validation."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def resolve_path(root: Path, value: str | Path) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else (root / path).resolve()


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    root = project_root()
    selected = Path(path) if path else root / "config" / "ldl_cad.json"
    if not selected.is_absolute():
        selected = (Path.cwd() / selected).resolve()
    payload = json.loads(selected.read_text(encoding="utf-8"))
    required = {
        "project",
        "candidate_study_id",
        "disease_study_id",
        "open_targets_endpoint",
        "l2g_candidate_threshold",
        "candidate_p_threshold",
        "minimum_instrument_f",
        "disease_prior_inclusion",
        "disease_slab_sd",
        "edge_prior_inclusion",
        "edge_slab_sd",
        "core_thresholds",
        "minimum_network_coverage",
        "source_locus_window_bp",
        "data_directory",
        "bayesv2g_directory",
        "results_directory",
    }
    missing = required.difference(payload)
    if missing:
        raise ValueError(f"Configuration is missing required fields: {sorted(missing)}")
    for key in ("disease_prior_inclusion", "edge_prior_inclusion", "minimum_network_coverage"):
        value = float(payload[key])
        if not 0 < value < 1:
            raise ValueError(f"{key} must lie strictly between zero and one")
    for key in ("disease_slab_sd", "edge_slab_sd", "minimum_instrument_f"):
        if float(payload[key]) <= 0:
            raise ValueError(f"{key} must be strictly positive")
    if not 0 <= float(payload["l2g_candidate_threshold"]) <= 1:
        raise ValueError("l2g_candidate_threshold must lie in [0, 1]")
    if not 0 < float(payload["candidate_p_threshold"]) <= 1:
        raise ValueError("candidate_p_threshold must lie in (0, 1]")
    if int(payload["source_locus_window_bp"]) <= 0:
        raise ValueError("source_locus_window_bp must be a positive integer")
    thresholds = [int(item) for item in payload["core_thresholds"]]
    if not thresholds or any(item < 1 for item in thresholds):
        raise ValueError("core_thresholds must contain positive integers")
    payload["core_thresholds"] = sorted(set(thresholds))
    if payload["core_thresholds"] != [2, 3, 5]:
        raise ValueError("This frozen report requires core_thresholds [2, 3, 5]")
    if payload["candidate_study_id"] != "GCST90019512":
        raise ValueError("This frozen acquisition workflow supports candidate_study_id GCST90019512")
    if payload["disease_study_id"] != "GCST90132314":
        raise ValueError("This frozen acquisition workflow supports disease_study_id GCST90132314")
    payload["config_path"] = str(selected)
    return payload
