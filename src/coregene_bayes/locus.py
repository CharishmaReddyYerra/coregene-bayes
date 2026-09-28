"""Real CAD/liver SORT1 locus audit."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import norm

from coregene_bayes.coloc import coloc_abf
from coregene_bayes.mr import harmonize_alleles, wald_ratio


GENE_SYMBOLS = {
    "ENSG00000134243": "SORT1",
    "ENSG00000143126": "CELSR2",
    "ENSG00000134222": "PSRC1",
}


def _snv_key(position: pd.Series, first: pd.Series, second: pd.Series) -> pd.Series:
    one = first.astype(str).str.upper()
    two = second.astype(str).str.upper()
    low = one.where(one <= two, two)
    high = two.where(one <= two, one)
    return position.astype(int).astype(str) + ":" + low + "/" + high


def analyse_sort1_locus(liver_path: str, cad_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    liver = pd.read_csv(liver_path, sep="\t")
    cad = pd.read_csv(cad_path, sep="\t", dtype={"p_value": str})
    liver = liver.loc[
        liver["type"].eq("SNP")
        & liver["ref"].astype(str).str.fullmatch("[ACGTacgt]")
        & liver["alt"].astype(str).str.fullmatch("[ACGTacgt]")
    ].copy()
    cad = cad.loc[
        cad["effect_allele"].astype(str).str.fullmatch("[ACGTacgt]")
        & cad["other_allele"].astype(str).str.fullmatch("[ACGTacgt]")
    ].copy()
    liver["variant_key"] = _snv_key(liver["position"], liver["ref"], liver["alt"])
    cad["variant_key"] = _snv_key(
        cad["base_pair_location"], cad["effect_allele"], cad["other_allele"]
    )
    liver = liver.drop_duplicates(["gene_id", "variant_key"])
    cad = cad.drop_duplicates(["variant_key"])

    summaries: list[dict] = []
    details: list[pd.DataFrame] = []
    for gene_id, qtl in liver.groupby("gene_id", sort=True):
        aligned = qtl.merge(cad, on="variant_key", suffixes=("_qtl", "_cad"), validate="one_to_one")
        if len(aligned) < 2:
            continue
        coloc = coloc_abf(
            aligned["beta_qtl"].to_numpy(float),
            aligned["se"].to_numpy(float),
            aligned["beta_cad"].to_numpy(float),
            aligned["standard_error"].to_numpy(float),
        )
        prior_values = {}
        for variance in (0.01, 0.16):
            prior_values[variance] = coloc_abf(
                aligned["beta_qtl"].to_numpy(float),
                aligned["se"].to_numpy(float),
                aligned["beta_cad"].to_numpy(float),
                aligned["standard_error"].to_numpy(float),
                prior_variance_one=variance,
                prior_variance_two=variance,
            ).pp_h4
        shared_prior_values = {
            p12: coloc_abf(
                aligned["beta_qtl"].to_numpy(float),
                aligned["se"].to_numpy(float),
                aligned["beta_cad"].to_numpy(float),
                aligned["standard_error"].to_numpy(float),
                p12=p12,
            ).pp_h4
            for p12 in (1e-6, 1e-5, 1e-4)
        }

        harmonized_rows: list[tuple[int, float]] = []
        for index, row in aligned.iterrows():
            result = harmonize_alleles(
                row["alt"],
                row["ref"],
                row["effect_allele"],
                row["other_allele"],
                float(row["beta_cad"]),
            )
            if result.retained:
                harmonized_rows.append((index, result.aligned_outcome_beta))
        mr_beta = mr_se = mr_p = instrument_f = np.nan
        instrument = None
        if harmonized_rows:
            chosen_index, aligned_beta = max(
                harmonized_rows,
                key=lambda item: abs(float(aligned.loc[item[0], "beta_qtl"]) / float(aligned.loc[item[0], "se"])),
            )
            chosen = aligned.loc[chosen_index]
            estimate = wald_ratio(
                float(chosen["beta_qtl"]),
                float(chosen["se"]),
                aligned_beta,
                float(chosen["standard_error"]),
            )
            mr_beta = estimate.beta
            mr_se = estimate.se
            mr_p = 2 * norm.sf(abs(mr_beta / mr_se))
            instrument_f = estimate.f_statistic
            instrument = chosen.get("rsid_qtl") or chosen.get("rsid_cad") or chosen["variant_key"]

        top = aligned.iloc[coloc.top_shared_index]
        summaries.append(
            {
                "gene_id": gene_id,
                "gene_symbol": GENE_SYMBOLS.get(gene_id, gene_id),
                "shared_variant_count": len(aligned),
                "pp_h0": coloc.posterior["H0"],
                "pp_h1": coloc.posterior["H1"],
                "pp_h2": coloc.posterior["H2"],
                "pp_h3": coloc.posterior["H3"],
                "pp_h4": coloc.posterior["H4"],
                "pp_h4_given_h3_h4": coloc.posterior["H4"]
                / max(coloc.posterior["H3"] + coloc.posterior["H4"], np.finfo(float).tiny),
                "pp_h4_prior_variance_0_01": prior_values[0.01],
                "pp_h4_prior_variance_0_16": prior_values[0.16],
                "pp_h4_shared_prior_1e_6": shared_prior_values[1e-6],
                "pp_h4_shared_prior_1e_5": shared_prior_values[1e-5],
                "pp_h4_shared_prior_1e_4": shared_prior_values[1e-4],
                "top_shared_variant": top.get("rsid_qtl") or top.get("rsid_cad") or top["variant_key"],
                "top_shared_variant_probability_given_h4": coloc.shared_variant_probability[coloc.top_shared_index],
                "mr_instrument": instrument,
                "mr_beta": mr_beta,
                "mr_se": mr_se,
                "mr_p_value": mr_p,
                "instrument_f": instrument_f,
                "mr_method": "single strongest liver cis-eQTL Wald ratio",
                "coloc_method": "classic single-signal Wakefield ABF",
            }
        )
        detail = aligned[
            [
                "variant_key",
                "rsid_qtl",
                "beta_qtl",
                "se",
                "beta_cad",
                "standard_error",
            ]
        ].copy()
        detail.insert(0, "gene_symbol", GENE_SYMBOLS.get(gene_id, gene_id))
        detail.insert(0, "gene_id", gene_id)
        detail["shared_variant_probability_given_h4"] = coloc.shared_variant_probability
        details.append(detail)
    if not summaries or not details:
        raise ValueError(
            "No SORT1-region gene had at least two shared SNVs after CAD/liver filtering"
        )
    return pd.DataFrame(summaries), pd.concat(details, ignore_index=True)
