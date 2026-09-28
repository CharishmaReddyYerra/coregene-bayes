"""Auditable single-signal Wakefield-ABF colocalization for the locus audit."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import logsumexp


HYPOTHESES = ("H0", "H1", "H2", "H3", "H4")


@dataclass(frozen=True, slots=True)
class ColocResult:
    posterior: dict[str, float]
    top_shared_index: int
    shared_variant_probability: np.ndarray

    @property
    def pp_h4(self) -> float:
        return self.posterior["H4"]


def wakefield_log_abf(beta: np.ndarray, se: np.ndarray, prior_variance: float = 0.04) -> np.ndarray:
    beta_values = np.asarray(beta, dtype=float)
    se_values = np.asarray(se, dtype=float)
    if beta_values.ndim != 1 or beta_values.size == 0 or beta_values.shape != se_values.shape:
        raise ValueError("beta and se must be equal-length, non-empty vectors")
    if not np.all(np.isfinite(beta_values)) or not np.all(np.isfinite(se_values)) or np.any(se_values <= 0):
        raise ValueError("beta/se must be finite and standard errors positive")
    if not np.isfinite(prior_variance) or prior_variance <= 0:
        raise ValueError("prior_variance must be positive")
    variance = se_values**2
    shrinkage = prior_variance / (variance + prior_variance)
    return 0.5 * (np.log1p(-shrinkage) + shrinkage * (beta_values / se_values) ** 2)


def _log_distinct_sum(first: np.ndarray, second: np.ndarray) -> float:
    if first.size < 2:
        return float("-inf")
    prefix = np.full(first.size + 1, -np.inf)
    suffix = np.full(first.size + 1, -np.inf)
    for index in range(first.size):
        prefix[index + 1] = np.logaddexp(prefix[index], second[index])
    for index in range(first.size - 1, -1, -1):
        suffix[index] = np.logaddexp(suffix[index + 1], second[index])
    return float(logsumexp(first + np.logaddexp(prefix[:-1], suffix[1:])))


def coloc_abf(
    beta_one: np.ndarray,
    se_one: np.ndarray,
    beta_two: np.ndarray,
    se_two: np.ndarray,
    *,
    p1: float = 1e-4,
    p2: float = 1e-4,
    p12: float = 1e-5,
    prior_variance_one: float = 0.04,
    prior_variance_two: float = 0.04,
) -> ColocResult:
    """Classic H0--H4 calculation assuming at most one signal per trait."""

    if not (0 < p12 <= min(p1, p2) < 1):
        raise ValueError("Require 0 < p12 <= min(p1,p2) < 1")
    first = wakefield_log_abf(beta_one, se_one, prior_variance_one)
    second = wakefield_log_abf(beta_two, se_two, prior_variance_two)
    if first.shape != second.shape:
        raise ValueError("Both traits must describe the same aligned variants")
    shared = first + second
    log_weights = np.array(
        [
            0.0,
            np.log(p1) + logsumexp(first),
            np.log(p2) + logsumexp(second),
            np.log(p1) + np.log(p2) + _log_distinct_sum(first, second),
            np.log(p12) + logsumexp(shared),
        ]
    )
    probabilities = np.exp(log_weights - logsumexp(log_weights))
    shared_probability = np.exp(shared - logsumexp(shared))
    return ColocResult(
        posterior=dict(zip(HYPOTHESES, probabilities.astype(float), strict=True)),
        top_shared_index=int(np.argmax(shared_probability)),
        shared_variant_probability=shared_probability,
    )

