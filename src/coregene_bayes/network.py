"""Directed influence-graph construction and convergence summaries."""

from __future__ import annotations

from collections import defaultdict

import numpy as np
import pandas as pd


def assign_physical_source_loci(candidates: pd.DataFrame, window_bp: int = 1_000_000) -> pd.Series:
    """Cluster candidate genes into conservative physical source loci.

    The cluster is a proxy for a shared cis/LD region, not proof of LD
    independence. It prevents co-located source genes from being counted as
    independent convergent perturbations.
    """

    required = {"gene_id", "gene_chromosome", "gene_start", "gene_end"}
    missing = required.difference(candidates.columns)
    if missing:
        raise ValueError(f"Candidate table lacks source-locus fields: {sorted(missing)}")
    if window_bp < 0:
        raise ValueError("window_bp must be non-negative")
    frame = candidates.copy()
    frame["gene_start"] = pd.to_numeric(frame["gene_start"], errors="coerce")
    frame["gene_end"] = pd.to_numeric(frame["gene_end"], errors="coerce")
    if frame[["gene_chromosome", "gene_start", "gene_end"]].isna().any().any():
        raise ValueError("Candidate genomic locations cannot be missing")
    frame = frame.sort_values(["gene_chromosome", "gene_start", "gene_end", "gene_id"])

    mapping: dict[str, str] = {}
    for chromosome, group in frame.groupby("gene_chromosome", sort=True):
        cluster_number = 0
        cluster_end: int | None = None
        for row in group.itertuples(index=False):
            start = int(row.gene_start)
            end = int(row.gene_end)
            if cluster_end is None or start > cluster_end + window_bp:
                cluster_number += 1
                cluster_end = end
            else:
                cluster_end = max(cluster_end, end)
            mapping[str(row.gene_id)] = f"chr{chromosome}:cluster{cluster_number}"
    return candidates["gene_id"].astype(str).map(mapping)


def edge_qc_label(
    *,
    harmonization_status: str,
    instrument_f: float,
    minimum_f: float,
    edge_probability: float,
    fdr: float,
) -> str:
    if harmonization_status != "aligned":
        return "failed_harmonization"
    if not np.isfinite(instrument_f) or instrument_f < minimum_f:
        return "weak_instrument"
    if fdr <= 0.05 and edge_probability >= 0.5:
        return "supported_pleiotropy_ambiguous"
    if edge_probability >= 0.5:
        return "model_supported_not_fdr"
    return "indeterminate"


def network_coverage(
    edges: pd.DataFrame,
    source_locus_universe: set[str],
    target_gene_ids: pd.Series,
) -> pd.DataFrame:
    """Count tested source loci without treating unavailable pairs as nulls."""

    if not source_locus_universe:
        raise ValueError("source_locus_universe cannot be empty")
    retained = (
        edges["harmonization_retained"].fillna(False).astype(bool)
        if "harmonization_retained" in edges.columns
        else edges["harmonization_status"].eq("aligned")
    )
    tested = (
        edges.loc[retained]
        .groupby("target_gene_id")["source_locus_id"]
        .nunique()
    )
    result = pd.DataFrame({"gene_id": target_gene_ids.astype(str).drop_duplicates()})
    result["tested_source_loci"] = result["gene_id"].map(tested).fillna(0).astype(int)
    result["source_locus_universe"] = len(source_locus_universe)
    result["network_coverage"] = result["tested_source_loci"] / len(source_locus_universe)
    return result


def top_edges_for_plot(edges: pd.DataFrame, *, limit: int = 40) -> pd.DataFrame:
    retained = (
        edges["harmonization_retained"].fillna(False).astype(bool)
        if "harmonization_retained" in edges.columns
        else edges["harmonization_status"].eq("aligned")
    )
    eligible = edges.loc[
        retained & np.isfinite(edges["edge_probability"])
    ].copy()
    return eligible.sort_values(
        ["edge_probability", "mr_p_value"], ascending=[False, True]
    ).head(limit)
