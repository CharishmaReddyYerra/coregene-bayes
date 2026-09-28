"""Checksum-rich provenance utilities."""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy
import pandas
import scipy


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_record(path: str | Path, *, root: Path | None = None) -> dict[str, Any]:
    selected = Path(path).resolve()
    name = str(selected.relative_to(root)) if root and selected.is_relative_to(root) else str(selected)
    return {"path": name, "bytes": selected.stat().st_size, "sha256": sha256_file(selected)}


def write_json(payload: dict[str, Any], path: str | Path) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + ".part")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(target)


def runtime_record() -> dict[str, str]:
    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": numpy.__version__,
        "pandas": pandas.__version__,
        "scipy": scipy.__version__,
    }


def manifest(
    *,
    project: str,
    config: dict[str, Any],
    inputs: Iterable[Path],
    outputs: Iterable[Path],
    root: Path,
    notes: dict[str, Any] | None = None,
) -> dict[str, Any]:
    return {
        "project": project,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "runtime": runtime_record(),
        "config": config,
        "inputs": [file_record(path, root=root) for path in sorted(inputs)],
        "outputs": [file_record(path, root=root) for path in sorted(outputs)],
        "notes": notes or {},
    }


def write_sha256s(paths: Iterable[Path], destination: Path, *, relative_to: Path) -> None:
    rows = [f"{sha256_file(path)}  {path.relative_to(relative_to)}" for path in sorted(paths)]
    destination.write_text("\n".join(rows) + "\n", encoding="utf-8")

