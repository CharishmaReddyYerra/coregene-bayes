"""Bayesian evidence probabilities for directed gene-network convergence.

The functions here return probabilities under explicit statistical models and
conditional-independence assumptions.  They are deliberately described as
*evidence probabilities*: they are neither labels of biological truth nor the
probability that the MR assumptions are correct.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike
from scipy.special import expit, logit


def spike_slab_pip(
    beta: ArrayLike,
    se: ArrayLike,
    prior_inclusion: float = 0.10,
    slab_sd: float = 0.20,
) -> float | np.ndarray:
    """Return spike-and-slab posterior evidence-inclusion probabilities.

    The observation model is ``beta_hat | theta ~ Normal(theta, se**2)``.
    Under the spike, ``theta = 0``; under the slab,
    ``theta ~ Normal(0, slab_sd**2)``.  ``prior_inclusion`` is the prior mass
    on the slab.  The closed-form normal-normal Bayes factor is evaluated on
    the log scale for numerical stability.

    This probability quantifies support for a non-zero effect under the stated
    prior and likelihood.  It does not by itself establish a causal edge.
    Scalars produce a Python ``float``; broadcast array inputs produce an
    ``ndarray``.
    """

    estimates, standard_errors = _broadcast_effect_inputs(beta, se)
    prior = _validate_probability(prior_inclusion, "prior_inclusion")
    slab = _positive_scalar(slab_sd, "slab_sd")

    if prior == 0.0:
        posterior = np.zeros(estimates.shape, dtype=float)
    elif prior == 1.0:
        posterior = np.ones(estimates.shape, dtype=float)
    else:
        # BF_10 = N(beta; 0, se^2 + slab^2) / N(beta; 0, se^2).
        # Writing the quadratic term via z and a shrinkage fraction avoids
        # subtracting two potentially huge inverse variances.
        log_se2 = 2.0 * np.log(standard_errors)
        log_slab2 = 2.0 * np.log(slab)
        log_total_variance = np.logaddexp(log_se2, log_slab2)
        # Evaluate z^2 * slab_fraction as one log-domain term.  Computing the
        # factors separately can produce inf * 0 -> NaN for extreme yet finite
        # effects and very narrow slabs.
        with np.errstate(over="ignore", divide="ignore", under="ignore"):
            log_quadratic = (
                2.0 * np.log(np.abs(estimates))
                - log_se2
                + log_slab2
                - log_total_variance
            )
            quadratic = np.exp(log_quadratic)
            log_bayes_factor = 0.5 * (
                log_se2 - log_total_variance + quadratic
            )
        posterior = expit(logit(prior) + log_bayes_factor)

    if posterior.ndim == 0:
        return float(posterior)
    return np.asarray(posterior, dtype=float)


def prior_sensitivity(
    beta: ArrayLike,
    se: ArrayLike,
    *,
    prior_inclusions: Sequence[float] = (0.01, 0.05, 0.10, 0.25, 0.50),
    slab_sds: Sequence[float] = (0.10, 0.20, 0.50, 1.00),
    observation_ids: Sequence[object] | None = None,
) -> pd.DataFrame:
    """Evaluate evidence probabilities over a prior/slab sensitivity grid.

    The result is tidy: one row per observation and prior combination, with a
    column named ``evidence_inclusion_probability`` to avoid implying that the
    latent inclusion state is observed truth.
    """

    estimates, standard_errors = _broadcast_effect_inputs(beta, se)
    if estimates.ndim > 1:
        raise ValueError("beta and se must be scalar or one-dimensional")
    estimates = np.atleast_1d(estimates)
    standard_errors = np.atleast_1d(standard_errors)
    if estimates.size == 0:
        raise ValueError("beta and se must contain at least one observation")

    priors = _validated_grid(
        prior_inclusions, "prior_inclusions", probability=True
    )
    slabs = _validated_grid(slab_sds, "slab_sds", probability=False)
    if observation_ids is None:
        ids = list(range(estimates.size))
    else:
        ids = list(observation_ids)
        if len(ids) != estimates.size:
            raise ValueError(
                "observation_ids must have the same length as beta and se"
            )
        if pd.isna(pd.Series(ids, dtype=object)).any():
            raise ValueError("observation_ids must not contain missing values")

    records: list[dict[str, object]] = []
    for prior in priors:
        for slab in slabs:
            probabilities = np.atleast_1d(
                spike_slab_pip(
                    estimates,
                    standard_errors,
                    prior_inclusion=prior,
                    slab_sd=slab,
                )
            )
            records.extend(
                {
                    "observation_id": observation_id,
                    "beta": float(effect),
                    "se": float(error),
                    "prior_inclusion": prior,
                    "slab_sd": slab,
                    "evidence_inclusion_probability": float(probability),
                }
                for observation_id, effect, error, probability in zip(
                    ids,
                    estimates,
                    standard_errors,
                    probabilities,
                    strict=True,
                )
            )
    return pd.DataFrame.from_records(records)


def collapse_locus_probabilities(
    frame: pd.DataFrame,
    target_col: str = "target_gene_id",
    locus_col: str = "source_locus_id",
    prob_col: str = "edge_probability",
    *,
    output_col: str = "locus_probability",
) -> pd.DataFrame:
    """Conservatively collapse edge evidence to locus-level probabilities.

    For all edge-level probabilities ``q_ij`` from locus ``l`` into target
    gene ``j``, the returned probability is

    ``q_lj = max_i(q_ij)``.

    Multiple variants or source genes from one locus therefore contribute at
    most one Bernoulli opportunity in downstream convergence counts.  Taking
    the maximum deliberately avoids treating correlated, and sometimes
    identical, within-locus evidence as independent replication.  It is a
    conservative evidence-collapse rule rather than a joint posterior for the
    biological locus event.
    """

    if not isinstance(frame, pd.DataFrame):
        raise TypeError("frame must be a pandas DataFrame")
    for value, name in (
        (target_col, "target_col"),
        (locus_col, "locus_col"),
        (prob_col, "prob_col"),
        (output_col, "output_col"),
    ):
        if not isinstance(value, str) or not value:
            raise ValueError(f"{name} must be a non-empty string")
    if len({target_col, locus_col, prob_col}) != 3:
        raise ValueError("target_col, locus_col, and prob_col must be distinct")
    if output_col in (target_col, locus_col, "n_edges_collapsed"):
        raise ValueError(
            "output_col must differ from identifier and reserved count columns"
        )
    required = {target_col, locus_col, prob_col}
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"frame is missing required columns: {missing}")
    if frame[target_col].isna().any() or frame[locus_col].isna().any():
        raise ValueError("target and locus identifiers must not be missing")

    try:
        probabilities = pd.to_numeric(frame[prob_col], errors="raise").to_numpy(
            dtype=float
        )
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{prob_col} must be numeric") from exc
    _validate_probability_array(probabilities, prob_col)

    working = frame[[target_col, locus_col]].copy()
    working["__probability"] = probabilities
    records: list[dict[str, object]] = []
    for (target, locus), group in working.groupby(
        [target_col, locus_col], sort=False, observed=True
    ):
        values = group["__probability"].to_numpy(dtype=float)
        collapsed = float(np.max(values))
        records.append(
            {
                target_col: target,
                locus_col: locus,
                output_col: collapsed,
                "n_edges_collapsed": int(values.size),
                "collapse_rule": "maximum_edge_probability",
            }
        )

    columns = [target_col, locus_col, output_col, "n_edges_collapsed", "collapse_rule"]
    return pd.DataFrame.from_records(records, columns=columns)


def core_probability(
    disease_probability: float,
    locus_probabilities: ArrayLike,
    thresholds: Iterable[int] = (2, 3, 5),
    *,
    method: str = "exact",
    monte_carlo_samples: int = 200_000,
    random_seed: int = 0,
) -> dict[int, float]:
    """Compute ``P(D=1 and at least m incoming loci)`` for each threshold.

    Each locus is modeled as an independent Bernoulli variable with its given
    evidence probability, and disease evidence ``D`` is assumed independent of
    those variables conditional on the model.  ``method='exact'`` evaluates the
    Poisson-binomial distribution by dynamic programming.  The
    ``'monte_carlo'`` method uses a fixed, configurable seed and is therefore
    deterministic for identical arguments.  ``'auto'`` uses exact evaluation
    for at most 2,000 loci and Monte Carlo above that point.

    Returned values are model-based convergence evidence probabilities, not
    declarations that a gene truly is or is not a core gene.
    """

    disease = _validate_probability(disease_probability, "disease_probability")
    probabilities = _as_probability_vector(
        locus_probabilities, "locus_probabilities"
    )
    requested_thresholds = _validate_thresholds(thresholds)
    if not isinstance(method, str):
        raise TypeError("method must be a string")
    normalized_method = method.strip().lower().replace("-", "_")
    if normalized_method == "mc":
        normalized_method = "monte_carlo"
    if normalized_method not in {"exact", "monte_carlo", "auto"}:
        raise ValueError("method must be 'exact', 'monte_carlo', or 'auto'")
    if normalized_method == "auto":
        normalized_method = "exact" if probabilities.size <= 2_000 else "monte_carlo"

    if normalized_method == "exact":
        tails = _poisson_binomial_tails(probabilities)
    else:
        samples = _positive_integer(monte_carlo_samples, "monte_carlo_samples")
        seed = _nonnegative_integer(random_seed, "random_seed")
        tails = _monte_carlo_tails(probabilities, samples, seed)

    n_loci = probabilities.size
    output: dict[int, float] = {}
    for threshold in requested_thresholds:
        incoming_probability = 0.0 if threshold > n_loci else float(tails[threshold])
        output[threshold] = float(np.clip(disease * incoming_probability, 0.0, 1.0))
    return output


def poisson_binomial_tail_probabilities(
    locus_probabilities: ArrayLike,
) -> np.ndarray:
    """Return exact ``P(N >= m)`` for ``m=0, ..., number of loci``."""

    probabilities = _as_probability_vector(
        locus_probabilities, "locus_probabilities"
    )
    return _poisson_binomial_tails(probabilities)


# An explicit alias for callers that prefer the longer method name.
spike_slab_prior_sensitivity = prior_sensitivity


def _broadcast_effect_inputs(
    beta: ArrayLike, se: ArrayLike
) -> tuple[np.ndarray, np.ndarray]:
    try:
        estimates, standard_errors = np.broadcast_arrays(
            np.asarray(beta, dtype=float), np.asarray(se, dtype=float)
        )
    except (TypeError, ValueError) as exc:
        raise TypeError("beta and se must be numeric and broadcast-compatible") from exc
    if not np.all(np.isfinite(estimates)):
        raise ValueError("beta must contain only finite values")
    if not np.all(np.isfinite(standard_errors)):
        raise ValueError("se must contain only finite values")
    if np.any(standard_errors <= 0.0):
        raise ValueError("se must contain only strictly positive values")
    return estimates, standard_errors


def _validate_probability(value: float, name: str) -> float:
    result = _finite_scalar(value, name)
    if not 0.0 <= result <= 1.0:
        raise ValueError(f"{name} must lie in [0, 1]")
    return result


def _validate_probability_array(values: np.ndarray, name: str) -> None:
    if not np.all(np.isfinite(values)):
        raise ValueError(f"{name} must contain only finite values")
    if np.any((values < 0.0) | (values > 1.0)):
        raise ValueError(f"{name} values must lie in [0, 1]")


def _as_probability_vector(values: ArrayLike, name: str) -> np.ndarray:
    if isinstance(values, np.ndarray):
        candidate = values
    elif isinstance(values, (pd.Series, list, tuple)):
        candidate = values
    else:
        try:
            candidate = list(values)  # type: ignore[arg-type]
        except TypeError as exc:
            raise TypeError(f"{name} must be a one-dimensional iterable") from exc
    try:
        array = np.asarray(candidate, dtype=float)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{name} must be numeric") from exc
    if array.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional")
    _validate_probability_array(array, name)
    return array


def _validated_grid(
    values: Sequence[float], name: str, *, probability: bool
) -> list[float]:
    try:
        candidates = list(values)
    except TypeError as exc:
        raise TypeError(f"{name} must be an iterable") from exc
    if not candidates:
        raise ValueError(f"{name} must not be empty")
    validated = [
        _validate_probability(value, name)
        if probability
        else _positive_scalar(value, name)
        for value in candidates
    ]
    if len(set(validated)) != len(validated):
        raise ValueError(f"{name} must not contain duplicate values")
    return validated


def _validate_thresholds(thresholds: Iterable[int]) -> list[int]:
    try:
        candidates = list(thresholds)
    except TypeError as exc:
        raise TypeError("thresholds must be an iterable of integers") from exc
    if not candidates:
        raise ValueError("thresholds must not be empty")
    validated: list[int] = []
    for threshold in candidates:
        if isinstance(threshold, (bool, np.bool_)) or not isinstance(
            threshold, (int, np.integer)
        ):
            raise TypeError("each threshold must be an integer")
        value = int(threshold)
        if value < 0:
            raise ValueError("thresholds must be non-negative")
        if value in validated:
            raise ValueError("thresholds must not contain duplicates")
        validated.append(value)
    return validated


def _poisson_binomial_tails(probabilities: np.ndarray) -> np.ndarray:
    n_loci = probabilities.size
    mass = np.zeros(n_loci + 1, dtype=float)
    mass[0] = 1.0
    for index, probability in enumerate(probabilities):
        # Reverse traversal prevents updated cells from contaminating the next
        # Bernoulli recurrence.
        for count in range(index + 1, 0, -1):
            mass[count] = (
                mass[count] * (1.0 - probability)
                + mass[count - 1] * probability
            )
        mass[0] *= 1.0 - probability
    tails = np.cumsum(mass[::-1])[::-1]
    tails[0] = 1.0  # Remove accumulated floating-point drift at the identity.
    return np.clip(tails, 0.0, 1.0)


def _monte_carlo_tails(
    probabilities: np.ndarray, samples: int, random_seed: int
) -> np.ndarray:
    n_loci = probabilities.size
    if n_loci == 0:
        return np.ones(1, dtype=float)
    rng = np.random.default_rng(random_seed)
    histogram = np.zeros(n_loci + 1, dtype=np.int64)
    # Bound the temporary random matrix to roughly five million values.
    chunk_size = max(1, min(10_000, 5_000_000 // n_loci))
    remaining = samples
    while remaining:
        current = min(remaining, chunk_size)
        counts = np.sum(
            rng.random((current, n_loci)) < probabilities,
            axis=1,
            dtype=np.int64,
        )
        histogram += np.bincount(counts, minlength=n_loci + 1)
        remaining -= current
    mass = histogram.astype(float) / samples
    tails = np.cumsum(mass[::-1])[::-1]
    tails[0] = 1.0
    return np.clip(tails, 0.0, 1.0)


def _finite_scalar(value: float, name: str) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be a real number")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise TypeError(f"{name} must be a real number") from exc
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _positive_scalar(value: float, name: str) -> float:
    result = _finite_scalar(value, name)
    if result <= 0.0:
        raise ValueError(f"{name} must be strictly positive")
    return result


def _integer(value: int, name: str) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(
        value, (int, np.integer)
    ):
        raise TypeError(f"{name} must be an integer")
    return int(value)


def _positive_integer(value: int, name: str) -> int:
    result = _integer(value, name)
    if result <= 0:
        raise ValueError(f"{name} must be strictly positive")
    return result


def _nonnegative_integer(value: int, name: str) -> int:
    result = _integer(value, name)
    if result < 0:
        raise ValueError(f"{name} must be non-negative")
    return result


__all__ = [
    "collapse_locus_probabilities",
    "core_probability",
    "poisson_binomial_tail_probabilities",
    "prior_sensitivity",
    "spike_slab_pip",
    "spike_slab_prior_sensitivity",
]
