"""End-to-end real-data CoreGene-Bayes analysis."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from coregene_bayes import __version__
from coregene_bayes.config import load_config, project_root, resolve_path
from coregene_bayes.evaluation import panelapp_green_genes, validation_summary
from coregene_bayes.locus import analyse_sort1_locus
from coregene_bayes.model import (
    collapse_locus_probabilities,
    core_probability,
    spike_slab_pip,
)
from coregene_bayes.mr import bh_fdr, harmonize_alleles, wald_ratio
from coregene_bayes.network import (
    assign_physical_source_loci,
    edge_qc_label,
    network_coverage,
)
from coregene_bayes.provenance import manifest, sha256_file, write_json, write_sha256s
from coregene_bayes.reporting import (
    plot_convergence,
    plot_core_rankings,
    plot_directed_network,
    plot_sort1_coloc,
    write_html_report,
)


REQUIRED_INPUTS = {
    "candidates.csv",
    "candidate_locus_rows.csv",
    "eqtlgen_cis_instrument_candidates.tsv.gz",
    "eqtlgen_network_instruments.tsv.gz",
    "eqtlgen_trans_outcomes.tsv.gz",
    "cad_instrument_outcomes.tsv.gz",
    "cad_sort1_region.tsv.gz",
    "gtex_v8_liver_sort1_genes.tsv.gz",
    "cad_metadata.yaml",
    "panelapp_fh_v2_6.json",
    "source_manifest.json",
}


def _verify_snapshot(directory: Path) -> list[Path]:
    missing = sorted(name for name in REQUIRED_INPUTS if not (directory / name).exists())
    if missing:
        raise FileNotFoundError(
            f"Real-data snapshot is incomplete ({missing}); run scripts/fetch_real_data.py"
        )
    source_manifest = json.loads((directory / "source_manifest.json").read_text(encoding="utf-8"))
    records = source_manifest.get("files")
    if not isinstance(records, dict) or not records:
        raise ValueError("Source manifest has no checksum records")
    unhashed = (REQUIRED_INPUTS - {"source_manifest.json"}).difference(records)
    if unhashed:
        raise ValueError(f"Required snapshot files are absent from source manifest: {sorted(unhashed)}")
    for filename, metadata in records.items():
        path = directory / filename
        if not path.exists():
            raise FileNotFoundError(f"Manifest-listed source is missing: {filename}")
        expected = metadata.get("sha256") if isinstance(metadata, dict) else None
        if expected != sha256_file(path):
            raise ValueError(f"Source checksum mismatch: {filename}")
    return [directory / name for name in sorted(REQUIRED_INPUTS)]


def _approximate_standardized_effect(zscore: float, sample_size: float) -> tuple[float, float]:
    """Convert Z/N to an approximate standardized effect and SE.

    eQTLGen Phase-I public tables expose Z and per-pair N but not beta/SE.
    For inverse-normal expression this first-order representation uses
    beta=Z/sqrt(N), SE=1/sqrt(N). It is exact only under a standardized-design
    approximation and is explicitly labelled as such in every output.
    """

    z = float(zscore)
    n = float(sample_size)
    if not np.isfinite(z) or not np.isfinite(n) or n <= 2:
        raise ValueError("Z score must be finite and sample size greater than two")
    root_n = np.sqrt(n)
    return z / root_n, 1.0 / root_n


def _per_allele_to_genotype_sd(
    beta: float,
    standard_error: float,
    effect_allele_frequency: float,
) -> tuple[float, float, float]:
    """Put a per-allele GWAS association on a genotype-SD scale.

    The eQTLGen ``Z/sqrt(N)`` approximation is an association per standard
    deviation of genotype.  For a biallelic variant under Hardy-Weinberg
    equilibrium, the allele-dosage SD is ``sqrt(2p(1-p))``.  Scaling both the
    GWAS beta and SE by this value makes the Wald numerator and denominator
    use the same genotype unit.  Allele flipping does not change the variance.
    """

    estimate = float(beta)
    error = float(standard_error)
    frequency = float(effect_allele_frequency)
    if not np.isfinite(estimate):
        raise ValueError("GWAS beta must be finite")
    if not np.isfinite(error) or error <= 0:
        raise ValueError("GWAS standard error must be finite and positive")
    if not np.isfinite(frequency) or not 0 < frequency < 1:
        raise ValueError("effect-allele frequency must lie strictly between zero and one")
    genotype_sd = float(np.sqrt(2.0 * frequency * (1.0 - frequency)))
    return estimate * genotype_sd, error * genotype_sd, genotype_sd


def _deduplicate_gwas(frame: pd.DataFrame) -> pd.DataFrame:
    key = "rsid" if "rsid" in frame.columns else "variant_id"
    selected = frame.loc[frame[key].notna()].copy()
    duplicates = selected.duplicated(key, keep=False)
    if duplicates.any():
        columns = [column for column in ("beta", "standard_error", "effect_allele", "other_allele") if column in selected]
        conflicts = selected.loc[duplicates].groupby(key)[columns].nunique(dropna=False).gt(1).any(axis=1)
        if conflicts.any():
            selected = selected.loc[~selected[key].isin(conflicts[conflicts].index)].copy()
    return selected.drop_duplicates(key).rename(columns={key: "SNP"})


def analyse_gene_disease(
    candidates: pd.DataFrame,
    cis_path: Path,
    cad_path: Path,
    *,
    prior_inclusion: float,
    slab_sd: float,
    minimum_f: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cis = pd.read_csv(cis_path, sep="\t")
    cad = _deduplicate_gwas(pd.read_csv(cad_path, sep="\t", dtype={"p_value": str}))
    merged = cis.merge(cad, on="SNP", how="left", suffixes=("_eqtl", "_cad"), validate="many_to_one")
    if "effect_allele_frequency" not in merged.columns:
        raise ValueError(
            "CAD GWAS must provide effect_allele_frequency so per-allele effects can be "
            "matched to the genotype-standardized eQTL scale"
        )
    candidate_lookup = candidates.set_index("gene_id")
    audit_rows: list[dict[str, Any]] = []
    usable_rows: list[dict[str, Any]] = []
    for row in merged.itertuples(index=False):
        if pd.isna(row.effect_allele) or pd.isna(row.other_allele) or pd.isna(row.beta):
            audit_rows.append(
                {
                    "analysis": "gene_to_cad",
                    "source_gene_id": row.source_gene_id,
                    "source_gene_symbol": row.GeneSymbol,
                    "target": "CAD",
                    "instrument_snp": row.SNP,
                    "exposure_effect_allele": row.AssessedAllele,
                    "exposure_other_allele": row.OtherAllele,
                    "outcome_effect_allele": np.nan,
                    "outcome_other_allele": np.nan,
                    "harmonization_status": "missing_outcome",
                    "harmonization_retained": False,
                    "analysis_retained": False,
                    "effect_scale_status": "not_evaluated",
                }
            )
            continue
        result = harmonize_alleles(
            row.AssessedAllele,
            row.OtherAllele,
            row.effect_allele,
            row.other_allele,
            float(row.beta),
        )
        record = {
            "analysis": "gene_to_cad",
            "source_gene_id": row.source_gene_id,
            "source_gene_symbol": row.GeneSymbol,
            "target": "CAD",
            "instrument_snp": row.SNP,
            "exposure_effect_allele": row.AssessedAllele,
            "exposure_other_allele": row.OtherAllele,
            "outcome_effect_allele": row.effect_allele,
            "outcome_other_allele": row.other_allele,
            "harmonization_status": result.status,
            "harmonization_retained": result.retained,
            "analysis_retained": False,
            "effect_scale_status": "not_evaluated",
        }
        audit_rows.append(record)
        if not result.retained:
            continue
        try:
            beta_y_scaled, se_y_scaled, genotype_sd = _per_allele_to_genotype_sd(
                result.aligned_beta,
                float(row.standard_error),
                float(row.effect_allele_frequency),
            )
        except (TypeError, ValueError):
            record["effect_scale_status"] = "missing_or_invalid_outcome_eaf"
            continue
        record["analysis_retained"] = True
        record["effect_scale_status"] = "per_allele_to_genotype_sd"
        record["outcome_effect_allele_frequency"] = float(row.effect_allele_frequency)
        beta_x, se_x = _approximate_standardized_effect(row.Zscore, row.NrSamples)
        estimate = wald_ratio(beta_x, se_x, beta_y_scaled, se_y_scaled)
        usable_rows.append(
            {
                **record,
                "gene_id": row.source_gene_id,
                "gene_symbol": row.GeneSymbol,
                "beta_x_approx": beta_x,
                "se_x_approx": se_x,
                "beta_y_aligned_per_allele": result.aligned_beta,
                "se_y_per_allele": float(row.standard_error),
                "beta_y_per_genotype_sd": beta_y_scaled,
                "se_y_per_genotype_sd": se_y_scaled,
                "genotype_sd": genotype_sd,
                "mr_beta": estimate.beta,
                "mr_se": estimate.se,
                "mr_ci_lower": estimate.ci_lower,
                "mr_ci_upper": estimate.ci_upper,
                "mr_p_value": estimate.p_value,
                "instrument_f": estimate.f_statistic,
                "eqtl_sample_size": int(row.NrSamples),
                "gwas_sample_size": row.sample_size if hasattr(row, "sample_size") else np.nan,
            }
        )

    usable = pd.DataFrame(usable_rows)
    if usable.empty:
        raise ValueError("No eQTL instruments could be aligned to the CAD GWAS")
    selected = (
        usable.sort_values(["gene_id", "instrument_f", "mr_p_value"], ascending=[True, False, True])
        .drop_duplicates("gene_id")
        .copy()
    )
    selected["mr_fdr"] = bh_fdr(selected["mr_p_value"].to_numpy(float))
    selected["disease_probability"] = spike_slab_pip(
        selected["mr_beta"].to_numpy(float),
        selected["mr_se"].to_numpy(float),
        prior_inclusion=prior_inclusion,
        slab_sd=slab_sd,
    )
    selected["qc_status"] = np.where(
        selected["instrument_f"] >= minimum_f,
        "pass_single_instrument_screen",
        "weak_instrument",
    )
    selected.loc[selected["instrument_f"] < minimum_f, "disease_probability"] = np.nan
    selected["mr_method"] = "single strongest available cis-eQTL Wald ratio (full delta SE)"
    selected["effect_scale"] = (
        "eQTL beta=Z/sqrt(N); CAD per-allele beta/SE multiplied by sqrt(2p(1-p))"
    )

    missing_rows = []
    observed = set(selected["gene_id"].astype(str))
    for gene_id, row in candidate_lookup.iterrows():
        if str(gene_id) in observed:
            continue
        missing_rows.append(
            {
                "gene_id": gene_id,
                "gene_symbol": row["gene_symbol"],
                "qc_status": "no_harmonizable_cis_instrument_in_cad",
                "disease_probability": np.nan,
                "mr_method": "not estimable",
                "effect_scale": "not estimable",
            }
        )
    if missing_rows:
        selected = pd.concat([selected, pd.DataFrame(missing_rows)], ignore_index=True, sort=False)
    return selected.sort_values(["qc_status", "gene_symbol"]).reset_index(drop=True), pd.DataFrame(audit_rows)


def analyse_gene_network(
    candidates: pd.DataFrame,
    instrument_path: Path,
    outcome_path: Path,
    *,
    prior_inclusion: float,
    slab_sd: float,
    minimum_f: float,
    source_locus_window_bp: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, set[str]]:
    instruments = pd.read_csv(instrument_path, sep="\t")
    outcomes = pd.read_csv(outcome_path, sep="\t")
    locus_mapping = dict(
        zip(
            candidates["gene_id"].astype(str),
            assign_physical_source_loci(candidates, source_locus_window_bp).astype(str),
            strict=True,
        )
    )
    candidate_symbol = candidates.set_index("gene_id")["gene_symbol"].astype(str).to_dict()
    merged = instruments.merge(outcomes, on="SNP", suffixes=("_source", "_target"), validate="many_to_many")
    merged = merged.loc[merged["source_gene_id"].astype(str) != merged["target_gene_id"].astype(str)].copy()
    records: list[dict[str, Any]] = []
    audit: list[dict[str, Any]] = []
    for row in merged.itertuples(index=False):
        result = harmonize_alleles(
            row.AssessedAllele_source,
            row.OtherAllele_source,
            row.AssessedAllele_target,
            row.OtherAllele_target,
            _approximate_standardized_effect(row.Zscore_target, row.NrSamples_target)[0],
        )
        source_id = str(row.source_gene_id)
        target_id = str(row.target_gene_id)
        audit_record = {
            "analysis": "gene_to_gene",
            "source_gene_id": source_id,
            "source_gene_symbol": row.GeneSymbol_source,
            "target_gene_id": target_id,
            "target_gene_symbol": row.GeneSymbol_target,
            "instrument_snp": row.SNP,
            "harmonization_status": result.status,
            "harmonization_retained": result.retained,
        }
        audit.append(audit_record)
        if not result.retained:
            records.append(audit_record)
            continue
        beta_x, se_x = _approximate_standardized_effect(row.Zscore_source, row.NrSamples_source)
        _, se_y = _approximate_standardized_effect(row.Zscore_target, row.NrSamples_target)
        estimate = wald_ratio(beta_x, se_x, result.aligned_beta, se_y)
        records.append(
            {
                **audit_record,
                "source_locus_id": locus_mapping[source_id],
                "beta_x_approx": beta_x,
                "se_x_approx": se_x,
                "beta_y_aligned_approx": result.aligned_beta,
                "se_y_approx": se_y,
                "mr_beta": estimate.beta,
                "mr_se": estimate.se,
                "mr_ci_lower": estimate.ci_lower,
                "mr_ci_upper": estimate.ci_upper,
                "mr_p_value": estimate.p_value,
                "instrument_f": estimate.f_statistic,
                "source_sample_size": int(row.NrSamples_source),
                "target_sample_size": int(row.NrSamples_target),
            }
        )
    edges = pd.DataFrame(records)
    if edges.empty:
        raise ValueError("No candidate gene-to-gene pairs were available after joining instruments and outcomes")
    valid = edges["harmonization_retained"].fillna(False) & edges["mr_p_value"].notna()
    edges["mr_fdr"] = np.nan
    edges.loc[valid, "mr_fdr"] = bh_fdr(edges.loc[valid, "mr_p_value"].to_numpy(float))
    edges["edge_probability"] = np.nan
    strong = valid & edges["instrument_f"].ge(minimum_f)
    edges.loc[strong, "edge_probability"] = spike_slab_pip(
        edges.loc[strong, "mr_beta"].to_numpy(float),
        edges.loc[strong, "mr_se"].to_numpy(float),
        prior_inclusion=prior_inclusion,
        slab_sd=slab_sd,
    )
    edges["edge_qc"] = [
        edge_qc_label(
            harmonization_status="aligned" if bool(retained) else str(status),
            instrument_f=float(instrument_f) if pd.notna(instrument_f) else np.nan,
            minimum_f=minimum_f,
            edge_probability=float(probability) if pd.notna(probability) else 0.0,
            fdr=float(fdr) if pd.notna(fdr) else 1.0,
        )
        for retained, status, instrument_f, probability, fdr in zip(
            edges["harmonization_retained"].fillna(False),
            edges["harmonization_status"],
            edges.get("instrument_f", pd.Series(np.nan, index=edges.index)),
            edges["edge_probability"],
            edges["mr_fdr"],
            strict=True,
        )
    ]
    edges["mr_method"] = "single cis instrument Wald ratio; same-cohort source/target summary data"
    edges["interpretation"] = "genetically anchored directed influence; horizontal pleiotropy unresolved"
    source_loci = {locus_mapping[str(gene)] for gene in instruments["source_gene_id"].astype(str) if str(gene) in locus_mapping}
    coverage = network_coverage(edges, source_loci, candidates["gene_id"])
    collapsed = collapse_locus_probabilities(
        edges.loc[strong, ["target_gene_id", "source_locus_id", "edge_probability"]],
        target_col="target_gene_id",
        locus_col="source_locus_id",
        prob_col="edge_probability",
    )
    return edges, collapsed, coverage, source_loci


def _build_core_table(
    candidates: pd.DataFrame,
    disease: pd.DataFrame,
    edges: pd.DataFrame,
    collapsed: pd.DataFrame,
    coverage: pd.DataFrame,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    disease_map = disease.set_index("gene_id")["disease_probability"].to_dict()
    collapsed_groups = {
        str(gene): group["locus_probability"].to_numpy(float)
        for gene, group in collapsed.groupby("target_gene_id")
    }
    rows: list[dict[str, Any]] = []
    thresholds = tuple(config["core_thresholds"])
    coverage_map = coverage.set_index("gene_id").to_dict("index")
    for candidate in candidates.itertuples(index=False):
        gene_id = str(candidate.gene_id)
        disease_probability = disease_map.get(gene_id, np.nan)
        locus_probabilities = collapsed_groups.get(gene_id, np.array([], dtype=float))
        core = (
            core_probability(float(disease_probability), locus_probabilities, thresholds)
            if np.isfinite(disease_probability)
            else {threshold: np.nan for threshold in thresholds}
        )
        coverage_record = coverage_map.get(gene_id, {})
        network_coverage_value = float(coverage_record.get("network_coverage", 0.0))
        record: dict[str, Any] = {
            "gene_id": gene_id,
            "gene_symbol": candidate.gene_symbol,
            "l2g_score_for_ascertainment": float(candidate.l2g_score),
            "bayesv2g_seed": bool(candidate.bayesv2g_seed),
            "disease_probability": disease_probability,
            "expected_incoming_loci": float(locus_probabilities.sum()),
            "tested_source_loci": int(coverage_record.get("tested_source_loci", 0)),
            "source_locus_universe": int(coverage_record.get("source_locus_universe", 0)),
            "network_coverage": network_coverage_value,
            "core_probability_scope": "conditional on observed/tested eQTLGen source loci and configured priors",
        }
        for threshold in thresholds:
            record[f"p_core_m{threshold}_observed"] = core[threshold]
            record[f"p_core_m{threshold}"] = (
                core[threshold]
                if network_coverage_value >= float(config["minimum_network_coverage"])
                else np.nan
            )
        rows.append(record)
    core = pd.DataFrame(rows)

    sensitivity_records: list[dict[str, Any]] = []
    disease_valid = disease.loc[
        disease["disease_probability"].notna()
        & disease["mr_beta"].notna()
        & disease["mr_se"].notna()
    ].copy()
    edge_valid = edges.loc[edges["edge_probability"].notna()].copy()
    disease_priors = sorted({0.05, float(config["disease_prior_inclusion"]), 0.20})
    disease_slabs = sorted({0.10, float(config["disease_slab_sd"]), 0.40})
    edge_priors = sorted({0.001, float(config["edge_prior_inclusion"]), 0.05})
    edge_slabs = sorted({0.075, float(config["edge_slab_sd"]), 0.30})
    edge_scenarios: list[tuple[float, float, dict[str, np.ndarray]]] = []
    for edge_prior in edge_priors:
        for edge_slab in edge_slabs:
            edge_valid["scenario_probability"] = spike_slab_pip(
                edge_valid["mr_beta"].to_numpy(float),
                edge_valid["mr_se"].to_numpy(float),
                prior_inclusion=edge_prior,
                slab_sd=edge_slab,
            )
            scenario_collapsed = collapse_locus_probabilities(
                edge_valid[["target_gene_id", "source_locus_id", "scenario_probability"]],
                prob_col="scenario_probability",
            )
            scenario_groups = {
                str(gene): group["locus_probability"].to_numpy(float)
                for gene, group in scenario_collapsed.groupby("target_gene_id")
            }
            edge_scenarios.append((edge_prior, edge_slab, scenario_groups))
    for disease_prior in disease_priors:
        for disease_slab in disease_slabs:
            disease_valid["scenario_probability"] = spike_slab_pip(
                disease_valid["mr_beta"].to_numpy(float),
                disease_valid["mr_se"].to_numpy(float),
                prior_inclusion=disease_prior,
                slab_sd=disease_slab,
            )
            scenario_disease = disease_valid.set_index("gene_id")["scenario_probability"].to_dict()
            for edge_prior, edge_slab, scenario_groups in edge_scenarios:
                for gene_id in core["gene_id"]:
                    disease_probability = scenario_disease.get(gene_id, np.nan)
                    probabilities = scenario_groups.get(gene_id, np.array([], dtype=float))
                    values = (
                        core_probability(float(disease_probability), probabilities, thresholds)
                        if np.isfinite(disease_probability)
                        else {threshold: np.nan for threshold in thresholds}
                    )
                    for threshold in thresholds:
                        sensitivity_records.append(
                            {
                                "gene_id": gene_id,
                                "disease_prior_inclusion": disease_prior,
                                "disease_slab_sd": disease_slab,
                                "edge_prior_inclusion": edge_prior,
                                "edge_slab_sd": edge_slab,
                                "threshold_m": threshold,
                                "core_evidence_probability": values[threshold],
                            }
                        )
    sensitivity = pd.DataFrame(sensitivity_records)
    primary = sensitivity.loc[sensitivity["threshold_m"].eq(2)].groupby("gene_id")[
        "core_evidence_probability"
    ].agg(["min", "max"])
    core["p_core_m2_sensitivity_min"] = core["gene_id"].map(primary["min"])
    core["p_core_m2_sensitivity_max"] = core["gene_id"].map(primary["max"])
    core = core.sort_values("p_core_m2", ascending=False, na_position="last").reset_index(drop=True)
    core["core_rank"] = np.where(core["p_core_m2"].notna(), core["p_core_m2"].rank(ascending=False, method="first"), np.nan)
    return core, sensitivity


def run_pipeline(
    *,
    config_path: str | Path | None = None,
    data_directory: str | Path | None = None,
    output_directory: str | Path | None = None,
) -> dict[str, Path]:
    root = project_root()
    config = load_config(config_path)
    data = Path(data_directory).resolve() if data_directory else resolve_path(root, config["data_directory"])
    output = Path(output_directory).resolve() if output_directory else resolve_path(root, config["results_directory"])
    output.mkdir(parents=True, exist_ok=True)
    figures = output / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    inputs = _verify_snapshot(data)

    candidates = pd.read_csv(data / "candidates.csv")
    disease, disease_audit = analyse_gene_disease(
        candidates,
        data / "eqtlgen_cis_instrument_candidates.tsv.gz",
        data / "cad_instrument_outcomes.tsv.gz",
        prior_inclusion=float(config["disease_prior_inclusion"]),
        slab_sd=float(config["disease_slab_sd"]),
        minimum_f=float(config["minimum_instrument_f"]),
    )
    edges, collapsed, coverage, source_loci = analyse_gene_network(
        candidates,
        data / "eqtlgen_network_instruments.tsv.gz",
        data / "eqtlgen_trans_outcomes.tsv.gz",
        prior_inclusion=float(config["edge_prior_inclusion"]),
        slab_sd=float(config["edge_slab_sd"]),
        minimum_f=float(config["minimum_instrument_f"]),
        source_locus_window_bp=int(config["source_locus_window_bp"]),
    )
    core, sensitivity = _build_core_table(candidates, disease, edges, collapsed, coverage, config)
    sort1, sort1_variants = analyse_sort1_locus(
        str(data / "gtex_v8_liver_sort1_genes.tsv.gz"),
        str(data / "cad_sort1_region.tsv.gz"),
    )
    validation = panelapp_green_genes(data / "panelapp_fh_v2_6.json")
    ranked, validation_metrics = validation_summary(core, validation)
    attrition = pd.DataFrame(
        [
            {
                "stage": "prespecified LDL-locus candidate genes",
                "count": len(candidates),
                "note": "Open Targets L2G rule plus verified BayesV2G seeds",
            },
            {
                "stage": "genes with a harmonized CAD Wald ratio",
                "count": int(disease["mr_beta"].notna().sum()),
                "note": "separate CAD GWAS; valid allele-frequency scaling required",
            },
            {
                "stage": "genes passing disease instrument/effect QC",
                "count": int(disease["disease_probability"].notna().sum()),
                "note": "F >= configured threshold",
            },
            {
                "stage": "candidate network source genes",
                "count": int(edges["source_gene_id"].nunique()),
                "note": "cis instrument present in fixed eQTLGen trans-tested SNP universe",
            },
            {
                "stage": "candidate network target genes",
                "count": int(edges["target_gene_id"].nunique()),
                "note": "complete tested trans coefficient available",
            },
            {
                "stage": "targets meeting network coverage threshold",
                "count": int(
                    coverage["network_coverage"]
                    .ge(float(config["minimum_network_coverage"]))
                    .sum()
                ),
                "note": "missing pairs were not coded as null edges",
            },
            {
                "stage": "genes eligible for primary core ranking",
                "count": int(ranked["p_core_m2"].notna().sum()),
                "note": "requires disease evidence and sufficient network coverage",
            },
        ]
    )

    tables = {
        "candidate_genes.csv": candidates,
        "gene_disease_mr.csv": disease,
        "gene_gene_edges.csv.gz": edges,
        "locus_collapsed_edges.csv": collapsed,
        "network_coverage.csv": coverage,
        "gene_features.csv": core.drop(columns=[column for column in core if column.startswith("p_core_")]),
        "core_gene_posteriors.csv": ranked,
        "prior_sensitivity.csv.gz": sensitivity,
        "sort1_liver_cad_colocalization.csv": sort1,
        "sort1_liver_cad_shared_variants.csv.gz": sort1_variants,
        "harmonization_audit.csv.gz": disease_audit,
        "validation_panelapp_genes.csv": validation,
        "validation_metrics.csv": validation_metrics,
        "analysis_attrition.csv": attrition,
    }
    written: dict[str, Path] = {}
    for filename, frame in tables.items():
        path = output / filename
        frame.to_csv(path, index=False, compression="gzip" if filename.endswith(".gz") else None)
        written[filename] = path

    plot_core_rankings(ranked, figures / "core_rankings.png")
    plot_convergence(ranked, figures / "convergence.png")
    plot_directed_network(edges, figures / "directed_network.png")
    plot_sort1_coloc(sort1, figures / "sort1_coloc.png")
    for path in figures.iterdir():
        if path.is_file():
            written[str(path.relative_to(output))] = path

    prior_locus_probabilities = np.full(len(source_loci), float(config["edge_prior_inclusion"]))
    prior_only_convergence = core_probability(
        1.0, prior_locus_probabilities, tuple(config["core_thresholds"])
    )
    prior_only_core = core_probability(
        float(config["disease_prior_inclusion"]),
        prior_locus_probabilities,
        tuple(config["core_thresholds"]),
    )
    model_payload = {
        "model_type": "spike-and-slab latent ratio-evidence-and-convergence event model",
        "version": __version__,
        "definition": "C_j(m) = D_j * 1(sum_l E_lj >= m)",
        "disease_prior_inclusion": config["disease_prior_inclusion"],
        "disease_slab_sd": config["disease_slab_sd"],
        "edge_prior_inclusion": config["edge_prior_inclusion"],
        "edge_slab_sd": config["edge_slab_sd"],
        "within_locus_collapse_rule": "maximum edge probability (conservative correlated-evidence rule)",
        "prior_only_convergence_probability": {
            f"m{threshold}": value for threshold, value in prior_only_convergence.items()
        },
        "prior_only_core_probability": {
            f"m{threshold}": value for threshold, value in prior_only_core.items()
        },
        "thresholds": config["core_thresholds"],
        "conditional_independence_assumptions": [
            "disease and incoming-edge latent indicators are independent conditional on model parameters",
            "collapsed physical source loci are treated as Bernoulli-independent",
        ],
    }
    write_json(model_payload, output / "model.json")
    written["model.json"] = output / "model.json"
    summary = {
        "candidate_gene_count": len(candidates),
        "gene_disease_estimate_count": int(disease["mr_beta"].notna().sum()),
        "directed_pair_count": len(edges),
        "source_locus_count": len(source_loci),
        "network_source_gene_count": int(edges["source_gene_id"].nunique()),
        "network_target_gene_count": int(edges["target_gene_id"].nunique()),
        "eligible_core_probability_count": int(ranked["p_core_m2"].notna().sum()),
        "top_gene": ranked.loc[ranked["p_core_m2"].idxmax(), "gene_symbol"] if ranked["p_core_m2"].notna().any() else None,
        "probability_scope": (
            "exploratory model-based evidence conditional on observed/tested pairs, candidate "
            "ascertainment, priors, and unverified instrumental-variable assumptions"
        ),
        "prior_only_convergence_probability_m2": prior_only_convergence[2],
        "prior_only_core_probability_m2": prior_only_core[2],
    }
    write_json(summary, output / "run_summary.json")
    written["run_summary.json"] = output / "run_summary.json"
    write_html_report(
        output / "report.html",
        core=ranked,
        edges=edges,
        disease_mr=disease,
        sort1=sort1,
        attrition=attrition,
        validation_metrics=validation_metrics,
        metadata={
            **config,
            **summary,
            "prior_only_convergence_probability": prior_only_convergence,
            "prior_only_core_probability": prior_only_core,
        },
    )
    written["report.html"] = output / "report.html"

    source_manifest = json.loads((data / "source_manifest.json").read_text(encoding="utf-8"))
    analysis_sources = [
        root / "pyproject.toml",
        root / "requirements.lock",
        Path(config["config_path"]),
        *sorted((root / "src" / "coregene_bayes").glob("*.py")),
    ]
    run_manifest = manifest(
        project=config["project"],
        config={key: value for key, value in config.items() if key != "config_path"},
        inputs=[*inputs, *analysis_sources],
        outputs=written.values(),
        root=root,
        notes={
            "source_snapshot_created_at": source_manifest.get("created_at_utc"),
            "source_counts": source_manifest.get("counts"),
            "summary": summary,
        },
    )
    write_json(run_manifest, output / "run_manifest.json")
    written["run_manifest.json"] = output / "run_manifest.json"
    write_sha256s(written.values(), output / "SHA256SUMS", relative_to=output)
    written["SHA256SUMS"] = output / "SHA256SUMS"
    return written
