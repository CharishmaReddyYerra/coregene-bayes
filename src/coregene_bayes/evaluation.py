"""External proxy-annotation overlap and ranking metrics."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


def panelapp_green_genes(path: str | Path) -> pd.DataFrame:
    payload = json.loads(Path(path).read_text(encoding="utf-8"))
    rows: list[dict[str, Any]] = []
    for gene in payload.get("genes", []):
        if str(gene.get("confidence_level")) != "3":
            continue
        data = gene.get("gene_data") or {}
        build = (data.get("ensembl_genes") or {}).get("GRch38", {})
        ensembl = next(iter(build.values()), {}).get("ensembl_id") if build else None
        rows.append(
            {
                "gene_symbol": data.get("gene_symbol") or gene.get("entity_name"),
                "gene_id": ensembl,
                "confidence_level": 3,
                "panel_id": payload.get("id"),
                "panel_name": payload.get("name"),
                "panel_version": payload.get("version"),
                "validation_scope": "diagnostic familial hypercholesterolaemia proxy, not core-gene truth",
            }
        )
    return pd.DataFrame(rows)


def average_precision(labels: np.ndarray, scores: np.ndarray) -> float:
    labels = np.asarray(labels, dtype=bool)
    scores = np.asarray(scores, dtype=float)
    valid = np.isfinite(scores)
    labels = labels[valid]
    scores = scores[valid]
    positives = int(labels.sum())
    if positives == 0:
        return float("nan")
    order = np.argsort(-scores, kind="stable")
    ranked = labels[order]
    precision = np.cumsum(ranked) / np.arange(1, ranked.size + 1)
    return float(precision[ranked].sum() / positives)


def validation_summary(
    ranking: pd.DataFrame,
    validation: pd.DataFrame,
    *,
    score_column: str = "p_core_m2",
    top_sizes: tuple[int, ...] = (10, 25, 50),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = ranking.copy()
    positive_ids = set(validation["gene_id"].dropna().astype(str))
    positive_symbols = set(validation["gene_symbol"].dropna().astype(str))
    frame["validation_positive"] = frame["gene_id"].astype(str).isin(positive_ids) | frame[
        "gene_symbol"
    ].astype(str).isin(positive_symbols)
    frame = frame.sort_values(score_column, ascending=False, na_position="last").reset_index(drop=True)
    frame["core_rank"] = np.where(
        frame[score_column].notna(),
        frame[score_column].rank(ascending=False, method="first"),
        np.nan,
    )
    eligible = frame[score_column].notna()
    universe = frame.loc[eligible]
    records: list[dict[str, Any]] = [
        {
            "metric": "panel_green_gene_count",
            "top_k": np.nan,
            "value": len(validation),
            "p_value": np.nan,
        },
        {
            "metric": "panel_candidate_overlap_count",
            "top_k": np.nan,
            "value": int(frame["validation_positive"].sum()),
            "p_value": np.nan,
        },
        {
            "metric": "panel_eligible_overlap_count",
            "top_k": np.nan,
            "value": int(universe["validation_positive"].sum()),
            "p_value": np.nan,
        },
        {
            "metric": "proxy_annotation_average_precision",
            "top_k": np.nan,
            "value": average_precision(
                universe["validation_positive"].to_numpy(), universe[score_column].to_numpy()
            ),
            "p_value": np.nan,
        }
    ]
    total_positive = int(universe["validation_positive"].sum())
    for top_k in top_sizes:
        selected = universe.head(min(top_k, len(universe)))
        a = int(selected["validation_positive"].sum())
        b = len(selected) - a
        c = total_positive - a
        d = len(universe) - len(selected) - c
        odds, pvalue = fisher_exact([[a, b], [c, d]], alternative="greater") if min(a, b, c, d) >= 0 else (np.nan, np.nan)
        records.append(
            {
                "metric": "proxy_annotation_top_k_enrichment",
                "top_k": top_k,
                "value": float(odds),
                "p_value": float(pvalue),
                "positive_in_top_k": a,
                "validation_positives_in_universe": total_positive,
                "eligible_gene_count": len(universe),
            }
        )
    return frame, pd.DataFrame(records)
