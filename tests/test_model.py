from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest
from scipy.stats import norm

from coregene_bayes.model import (
    collapse_locus_probabilities,
    core_probability,
    poisson_binomial_tail_probabilities,
    prior_sensitivity,
    spike_slab_pip,
)


def test_spike_slab_pip_matches_closed_form_bayes_factor() -> None:
    beta = 0.3
    se = 0.1
    prior = 0.2
    slab = 0.4
    null_density = norm.pdf(beta, loc=0.0, scale=se)
    slab_density = norm.pdf(beta, loc=0.0, scale=math.sqrt(se**2 + slab**2))
    bayes_factor = slab_density / null_density
    expected = prior * bayes_factor / (1.0 - prior + prior * bayes_factor)

    assert spike_slab_pip(beta, se, prior, slab) == pytest.approx(expected)


def test_spike_slab_pip_is_symmetric_vectorized_and_stable_for_large_z() -> None:
    result = spike_slab_pip(
        np.array([-0.5, 0.0, 0.5, 100.0]),
        0.1,
        prior_inclusion=0.1,
        slab_sd=0.3,
    )
    assert isinstance(result, np.ndarray)
    assert result.shape == (4,)
    assert result[0] == pytest.approx(result[2])
    assert result[1] < 0.1  # Occam penalty when the estimate is exactly null.
    assert result[3] == pytest.approx(1.0)
    assert np.all((result >= 0.0) & (result <= 1.0))
    extreme = spike_slab_pip(
        beta=1e60,
        se=1e-100,
        prior_inclusion=0.1,
        slab_sd=1e-300,
    )
    assert np.isfinite(extreme)
    assert 0.0 <= extreme <= 1.0


def test_spike_slab_prior_boundaries_and_validation() -> None:
    assert spike_slab_pip(1.0, 0.1, 0.0, 0.2) == 0.0
    assert spike_slab_pip(0.0, 0.1, 1.0, 0.2) == 1.0
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        spike_slab_pip(0.1, 0.1, prior_inclusion=1.1)
    with pytest.raises(ValueError, match="strictly positive"):
        spike_slab_pip(0.1, 0.0)
    with pytest.raises(ValueError, match="finite"):
        spike_slab_pip(np.nan, 0.1)


def test_prior_sensitivity_returns_tidy_evidence_probabilities() -> None:
    result = prior_sensitivity(
        beta=[0.1, 0.5],
        se=[0.1, 0.1],
        prior_inclusions=(0.05, 0.50),
        slab_sds=(0.2, 0.5),
        observation_ids=("edge-a", "edge-b"),
    )

    assert len(result) == 2 * 2 * 2
    assert set(result.columns) == {
        "observation_id",
        "beta",
        "se",
        "prior_inclusion",
        "slab_sd",
        "evidence_inclusion_probability",
    }
    assert result["evidence_inclusion_probability"].between(0.0, 1.0).all()
    paired = result[result["slab_sd"] == 0.2].pivot(
        index="observation_id",
        columns="prior_inclusion",
        values="evidence_inclusion_probability",
    )
    assert np.all(paired[0.50] > paired[0.05])
    assert "truth" not in " ".join(result.columns).lower()
    with pytest.raises(ValueError, match="at least one"):
        prior_sensitivity([], [])


def test_collapse_locus_probabilities_uses_conservative_maximum() -> None:
    frame = pd.DataFrame(
        {
            "source_gene_id": ["a", "b", "c", "d", "e"],
            "source_locus_id": ["l1", "l1", "l2", "l2", "l3"],
            "target_gene_id": ["target-1", "target-1", "target-1", "target-2", "target-2"],
            "edge_probability": [0.2, 0.3, 1.0, 0.0, 0.4],
        }
    )
    result = collapse_locus_probabilities(frame)

    assert result.columns.tolist() == [
        "target_gene_id",
        "source_locus_id",
        "locus_probability",
        "n_edges_collapsed",
        "collapse_rule",
    ]
    lookup = result.set_index(["target_gene_id", "source_locus_id"])
    assert lookup.loc[("target-1", "l1"), "locus_probability"] == pytest.approx(0.3)
    assert lookup.loc[("target-1", "l1"), "n_edges_collapsed"] == 2
    assert lookup.loc[("target-1", "l1"), "collapse_rule"] == "maximum_edge_probability"
    assert lookup.loc[("target-1", "l2"), "locus_probability"] == 1.0
    assert lookup.loc[("target-2", "l2"), "locus_probability"] == 0.0


def test_collapse_locus_probabilities_validates_input() -> None:
    with pytest.raises(ValueError, match="missing required"):
        collapse_locus_probabilities(pd.DataFrame({"edge_probability": [0.2]}))
    frame = pd.DataFrame(
        {
            "target_gene_id": ["g"],
            "source_locus_id": ["l"],
            "edge_probability": [1.2],
        }
    )
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        collapse_locus_probabilities(frame)
    with pytest.raises(ValueError, match="reserved"):
        collapse_locus_probabilities(frame, output_col="n_edges_collapsed")
    with pytest.raises(ValueError, match="distinct"):
        collapse_locus_probabilities(frame, target_col="source_locus_id")


def test_core_probability_exact_matches_hand_calculation() -> None:
    # For two fair Bernoulli loci, P(N >= 1)=.75 and P(N >= 2)=.25.
    result = core_probability(
        0.8, [0.5, 0.5], thresholds=(0, 1, 2, 3), method="exact"
    )
    assert result == pytest.approx({0: 0.8, 1: 0.6, 2: 0.2, 3: 0.0})


def test_exact_poisson_binomial_tails_cover_boundary_probabilities() -> None:
    tails = poisson_binomial_tail_probabilities([1.0, 0.0, 0.5])
    # One success is certain, and a second occurs with probability .5.
    np.testing.assert_allclose(tails, [1.0, 1.0, 0.5, 0.0])
    empty = poisson_binomial_tail_probabilities([])
    np.testing.assert_array_equal(empty, [1.0])


def test_core_probability_monte_carlo_is_deterministic_and_close_to_exact() -> None:
    probabilities = [0.15, 0.30, 0.55, 0.80]
    exact = core_probability(0.7, probabilities, thresholds=(1, 2, 3), method="exact")
    first = core_probability(
        0.7,
        probabilities,
        thresholds=(1, 2, 3),
        method="monte_carlo",
        monte_carlo_samples=60_000,
        random_seed=44,
    )
    second = core_probability(
        0.7,
        probabilities,
        thresholds=(1, 2, 3),
        method="monte_carlo",
        monte_carlo_samples=60_000,
        random_seed=44,
    )
    assert first == second
    for threshold in exact:
        assert first[threshold] == pytest.approx(exact[threshold], abs=0.006)


@pytest.mark.parametrize(
    ("arguments", "error", "message"),
    [
        ((-0.1, [0.2], (1,)), ValueError, "disease_probability"),
        ((0.5, [np.nan], (1,)), ValueError, "finite"),
        ((0.5, [1.2], (1,)), ValueError, r"\[0, 1\]"),
        ((0.5, [0.2], (-1,)), ValueError, "non-negative"),
        ((0.5, [0.2], (1.5,)), TypeError, "integer"),
        ((0.5, [0.2], (1, 1)), ValueError, "duplicates"),
    ],
)
def test_core_probability_validation(arguments, error, message) -> None:
    with pytest.raises(error, match=message):
        core_probability(*arguments)
