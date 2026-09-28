"""Validated file-level contract with the sibling BayesV2G project."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from coregene_bayes.provenance import sha256_file


REQUIRED_COLUMNS = {
    "studyLocusId",
    "geneId",
    "gene_symbol",
    "posterior_mean",
    "posterior_lower",
    "posterior_upper",
    "conditional_posterior_mean",
    "unresolved_probability",
    "rank",
    "probability_scope",
    "evidence_coverage",
}


def _parse_checksums(path: Path) -> dict[str, str]:
    checksums: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        digest, filename = line.split(maxsplit=1)
        checksums[filename.strip()] = digest
    return checksums


def load_bayesv2g_candidates(directory: str | Path) -> tuple[pd.DataFrame, dict]:
    """Load and checksum-verify the BayesV2G candidate artifact.

    BayesV2G probabilities define candidate ascertainment only. They are not
    reused as independent causal evidence in the CoreGene model because they
    already incorporate molecular-QTL colocalization features.
    """

    root = Path(directory).resolve()
    candidates_path = root / "prioritized_genes.csv"
    manifest_path = root / "run_manifest.json"
    checksums_path = root / "SHA256SUMS"
    for path in (candidates_path, manifest_path, checksums_path):
        if not path.exists():
            raise FileNotFoundError(f"Missing BayesV2G contract file: {path}")

    checksums = _parse_checksums(checksums_path)
    expected = checksums.get("prioritized_genes.csv")
    if expected is None:
        raise ValueError("BayesV2G SHA256SUMS does not cover prioritized_genes.csv")
    observed = sha256_file(candidates_path)
    if observed != expected:
        raise ValueError("BayesV2G candidate checksum does not match SHA256SUMS")

    candidates = pd.read_csv(candidates_path)
    missing = REQUIRED_COLUMNS.difference(candidates.columns)
    if missing:
        raise ValueError(f"BayesV2G candidates are missing columns: {sorted(missing)}")
    if candidates["geneId"].isna().any() or candidates["studyLocusId"].isna().any():
        raise ValueError("BayesV2G candidate identifiers cannot be missing")
    if candidates.duplicated(["studyLocusId", "geneId"]).any():
        raise ValueError("BayesV2G contains duplicate locus-gene rows")
    manifest_payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    return candidates, manifest_payload


def top_seed_genes(candidates: pd.DataFrame) -> pd.DataFrame:
    """Return one BayesV2G rank-one seed per configured locus."""

    seeds = candidates.loc[candidates["rank"].eq(1)].copy()
    return seeds.sort_values(["studyLocusId", "geneId"]).reset_index(drop=True)

