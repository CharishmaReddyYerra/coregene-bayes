from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from coregene_bayes.mr import (
    bh_fdr,
    harmonize_alleles,
    harmonize_dataframe,
    inverse_variance_weighted,
    mr_egger,
    wald_ratio,
    weighted_median,
)


@pytest.mark.parametrize(
    ("outcome_alleles", "expected_beta", "expected_status"),
    [
        (("A", "G"), 0.2, "aligned"),
        (("G", "A"), -0.2, "flipped"),
        (("T", "C"), 0.2, "strand_aligned"),
        (("C", "T"), -0.2, "strand_flipped"),
    ],
)
def test_harmonize_non_palindromic_signed_orientations(
    outcome_alleles, expected_beta, expected_status
) -> None:
    result = harmonize_alleles(
        "a", "g", outcome_alleles[0], outcome_alleles[1], 0.2
    )
    assert result.aligned_beta == pytest.approx(expected_beta)
    assert result.status == expected_status
    assert result.retained
    # The result intentionally supports tuple unpacking for pipeline use.
    aligned_beta, status = result
    assert aligned_beta == pytest.approx(expected_beta)
    assert status == expected_status


def test_harmonize_palindrome_uses_frequencies_or_drops_ambiguity() -> None:
    aligned = harmonize_alleles("A", "T", "A", "T", 0.3, 0.20, 0.21)
    flipped = harmonize_alleles("A", "T", "T", "A", 0.3, 0.20, 0.79)
    missing = harmonize_alleles("A", "T", "A", "T", 0.3)
    near_half = harmonize_alleles("C", "G", "C", "G", 0.3, 0.49, 0.51)

    assert aligned.aligned_beta == pytest.approx(0.3)
    assert aligned.status == "palindromic_aligned"
    assert flipped.aligned_beta == pytest.approx(-0.3)
    assert flipped.status == "palindromic_flipped"
    assert math.isnan(missing.aligned_beta)
    assert missing.status == "ambiguous_palindromic"
    assert not missing.retained
    assert math.isnan(near_half.aligned_beta)
    assert near_half.status == "ambiguous_palindromic"


def test_harmonize_incompatible_and_invalid_alleles() -> None:
    incompatible = harmonize_alleles("A", "G", "A", "C", 0.1)
    assert math.isnan(incompatible.aligned_beta)
    assert incompatible.status == "incompatible"

    with pytest.raises(ValueError, match="one of A, C, G, T"):
        harmonize_alleles("I", "D", "A", "G", 0.1)
    with pytest.raises(ValueError, match="must differ"):
        harmonize_alleles("A", "A", "A", "G", 0.1)
    with pytest.raises(ValueError, match="eaf_tolerance"):
        harmonize_alleles(
            "A", "T", "A", "T", 0.1, 0.2, 0.2, eaf_tolerance=0.5
        )


def test_harmonize_dataframe_preserves_audit_rows() -> None:
    frame = pd.DataFrame(
        {
            "exposure_effect_allele": ["A", "A", "A"],
            "exposure_other_allele": ["G", "G", "T"],
            "outcome_effect_allele": ["A", "G", "A"],
            "outcome_other_allele": ["G", "A", "T"],
            "beta_y": [0.2, 0.3, 0.4],
            "exposure_eaf": [0.2, 0.2, np.nan],
            "outcome_eaf": [0.2, 0.8, np.nan],
        }
    )
    result = harmonize_dataframe(frame)

    assert len(result) == len(frame)
    np.testing.assert_allclose(result.loc[:1, "beta_y_aligned"], [0.2, -0.3])
    assert np.isnan(result.loc[2, "beta_y_aligned"])
    assert result["harmonization_retained"].tolist() == [True, True, False]


def test_wald_ratio_uses_full_delta_standard_error_and_f_statistic() -> None:
    result = wald_ratio(beta_x=0.4, se_x=0.05, beta_y=0.2, se_y=0.04)

    expected_variance = 0.04**2 / 0.4**2 + 0.2**2 * 0.05**2 / 0.4**4
    assert result.beta == pytest.approx(0.5)
    assert result.se == pytest.approx(math.sqrt(expected_variance))
    assert result.f_statistic == pytest.approx(64.0)
    assert result.min_f_statistic == pytest.approx(64.0)
    assert result.ci_lower < result.beta < result.ci_upper
    assert 0.0 <= result.p_value <= 1.0


def test_wald_ratio_covariance_term_and_validation() -> None:
    independent = wald_ratio(0.4, 0.05, 0.2, 0.04)
    positively_correlated = wald_ratio(
        0.4, 0.05, 0.2, 0.04, covariance=0.001
    )
    assert positively_correlated.se < independent.se

    with pytest.raises(ValueError, match="non-zero"):
        wald_ratio(0.0, 0.05, 0.2, 0.04)
    with pytest.raises(ValueError, match="strictly positive"):
        wald_ratio(0.4, 0.0, 0.2, 0.04)
    with pytest.raises(ValueError, match="covariance"):
        wald_ratio(0.4, 0.05, 0.2, 0.04, covariance=0.003)


def test_ivw_recovers_common_ratio_and_reports_heterogeneity() -> None:
    result = inverse_variance_weighted(
        beta_x=[0.4, 0.2, -0.3],
        se_x=[0.05, 0.04, 0.05],
        beta_y=[0.2, 0.1, -0.15],
        se_y=[0.04, 0.03, 0.04],
    )

    assert result.beta == pytest.approx(0.5)
    assert result.q_statistic == pytest.approx(0.0, abs=1e-28)
    assert result.q_p_value == pytest.approx(1.0)
    assert result.n_instruments == 3
    assert result.se > 0.0


def test_ivw_random_effects_inflate_uncertainty_under_heterogeneity() -> None:
    arguments = dict(
        beta_x=[0.3, 0.3, 0.3],
        se_x=[0.02, 0.02, 0.02],
        beta_y=[0.03, 0.15, 0.30],
        se_y=[0.01, 0.01, 0.01],
    )
    fixed = inverse_variance_weighted(**arguments)
    random = inverse_variance_weighted(**arguments, random_effects=True)
    assert random.beta == pytest.approx(fixed.beta)
    assert random.se > fixed.se
    assert fixed.q_statistic > 0.0


def test_weighted_median_is_reproducible_and_resists_one_outlier() -> None:
    arguments = dict(
        beta_x=[1.0, 1.0, 1.0, 1.0, 1.0],
        se_x=[0.005] * 5,
        beta_y=[0.19, 0.20, 0.21, 0.22, 8.0],
        se_y=[0.05] * 5,
        n_bootstrap=1_500,
        random_seed=91,
    )
    first = weighted_median(**arguments)
    second = weighted_median(**arguments)

    assert first.beta == pytest.approx(0.21)
    assert first.se == second.se
    assert first.p_value == second.p_value
    assert first.beta < 0.3


def test_mr_egger_recovers_slope_and_intercept() -> None:
    exposure_increasing_beta = np.array([0.08, 0.14, 0.23, 0.31, 0.40])
    oriented_outcome_beta = 0.08 + 1.7 * exposure_increasing_beta
    # Arbitrarily reported effect alleles must not change Egger.  The method
    # re-orients both associations internally to beta_x > 0.
    reported_sign = np.array([-1.0, 1.0, -1.0, 1.0, 1.0])
    beta_x = exposure_increasing_beta * reported_sign
    beta_y = oriented_outcome_beta * reported_sign
    result = mr_egger(
        beta_x,
        np.full(beta_x.shape, 0.02),
        beta_y,
        np.full(beta_x.shape, 0.03),
        random_effects=False,
    )

    assert result.beta == pytest.approx(1.7, abs=1e-12)
    assert result.intercept == pytest.approx(0.08, abs=1e-12)
    assert result.intercept_se is not None and result.intercept_se > 0.0
    assert result.intercept_p_value is not None
    assert result.q_statistic == pytest.approx(0.0, abs=1e-26)


def test_multi_instrument_estimators_require_sufficient_data() -> None:
    with pytest.raises(ValueError, match="at least two"):
        inverse_variance_weighted([0.2], [0.02], [0.1], [0.03])
    with pytest.raises(ValueError, match="at least three"):
        weighted_median([0.2, 0.3], [0.02, 0.02], [0.1, 0.2], [0.03, 0.03])
    with pytest.raises(ValueError, match="at least three"):
        mr_egger([0.2, 0.3], [0.02, 0.02], [0.1, 0.2], [0.03, 0.03])
    with pytest.raises(ValueError, match="non-identical"):
        mr_egger(
            [0.2, 0.2, 0.2],
            [0.02] * 3,
            [0.1, 0.2, 0.3],
            [0.03] * 3,
        )


def test_bh_fdr_known_values_ties_and_missingness() -> None:
    values = np.array([0.01, 0.04, 0.03, 0.002, np.nan, 0.04])
    adjusted = bh_fdr(values)
    # Five observed hypotheses; the final p=.04 sets the reverse cumulative
    # minimum to .04 for the three largest ordered p-values.
    np.testing.assert_allclose(
        adjusted[[0, 1, 2, 3, 5]], [0.025, 0.04, 0.04, 0.01, 0.04]
    )
    assert np.isnan(adjusted[4])
    tied = bh_fdr([0.05, 0.05, 0.05])
    np.testing.assert_allclose(tied, 0.05)


@pytest.mark.parametrize("bad", [[-0.1], [1.01], [np.inf], [[0.1]]])
def test_bh_fdr_rejects_invalid_values(bad) -> None:
    with pytest.raises(ValueError):
        bh_fdr(bad)
