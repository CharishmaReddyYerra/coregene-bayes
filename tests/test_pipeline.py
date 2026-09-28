from __future__ import annotations

import math

import pandas as pd
import pytest

from coregene_bayes.pipeline import _per_allele_to_genotype_sd, analyse_gene_disease


def test_per_allele_to_genotype_sd_uses_real_allele_frequency() -> None:
    beta, se, genotype_sd = _per_allele_to_genotype_sd(0.20, 0.04, 0.25)
    expected_sd = math.sqrt(2.0 * 0.25 * 0.75)
    assert genotype_sd == pytest.approx(expected_sd)
    assert beta == pytest.approx(0.20 * expected_sd)
    assert se == pytest.approx(0.04 * expected_sd)


@pytest.mark.parametrize("frequency", [0.0, 1.0, float("nan")])
def test_per_allele_to_genotype_sd_rejects_invalid_frequency(frequency: float) -> None:
    with pytest.raises(ValueError, match="frequency"):
        _per_allele_to_genotype_sd(0.20, 0.04, frequency)


def test_gene_disease_uses_eaf_scale_and_audits_missing_outcome(tmp_path) -> None:
    candidates = pd.DataFrame({"gene_id": ["ENSG1"], "gene_symbol": ["GENE1"]})
    cis = pd.DataFrame(
        {
            "SNP": ["rs1", "rs_missing"],
            "AssessedAllele": ["A", "C"],
            "OtherAllele": ["G", "T"],
            "Zscore": [10.0, 8.0],
            "NrSamples": [100, 100],
            "GeneSymbol": ["GENE1", "GENE1"],
            "source_gene_id": ["ENSG1", "ENSG1"],
        }
    )
    cad = pd.DataFrame(
        {
            "rsid": ["rs1"],
            "effect_allele": ["A"],
            "other_allele": ["G"],
            "beta": [0.20],
            "standard_error": [0.04],
            "effect_allele_frequency": [0.25],
            "sample_size": [1000],
        }
    )
    cis_path = tmp_path / "cis.tsv.gz"
    cad_path = tmp_path / "cad.tsv.gz"
    cis.to_csv(cis_path, sep="\t", index=False, compression="gzip")
    cad.to_csv(cad_path, sep="\t", index=False, compression="gzip")

    estimates, audit = analyse_gene_disease(
        candidates,
        cis_path,
        cad_path,
        prior_inclusion=0.1,
        slab_sd=0.2,
        minimum_f=10.0,
    )

    expected_sd = math.sqrt(2 * 0.25 * 0.75)
    row = estimates.loc[estimates["gene_id"].eq("ENSG1")].iloc[0]
    assert row["beta_y_per_genotype_sd"] == pytest.approx(0.20 * expected_sd)
    assert row["se_y_per_genotype_sd"] == pytest.approx(0.04 * expected_sd)
    assert row["mr_beta"] == pytest.approx(0.20 * expected_sd)
    assert set(audit["harmonization_status"]) == {"aligned", "missing_outcome"}
