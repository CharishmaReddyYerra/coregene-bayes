"""Small, auditable Mendelian-randomization utilities.

The estimators in this module operate on *harmonized, approximately
independent* instruments.  They quantify evidence under the instrumental
variable assumptions; they do not establish biological truth.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import NamedTuple

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike
from scipy.stats import chi2, norm


_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
_PALINDROMIC = {frozenset(("A", "T")), frozenset(("C", "G"))}


class HarmonizationResult(NamedTuple):
    """Result of orienting an outcome effect to the exposure effect allele.

    ``aligned_beta`` is ``nan`` when ``status`` indicates that the variant
    cannot be retained.  A named tuple keeps the convenient ``beta, status =``
    unpacking API while also providing explicit attribute names.
    """

    aligned_beta: float
    status: str

    @property
    def retained(self) -> bool:
        return bool(np.isfinite(self.aligned_beta))

    @property
    def aligned_outcome_beta(self) -> float:
        """Backward-compatible explicit alias for :attr:`aligned_beta`."""

        return self.aligned_beta


@dataclass(frozen=True)
class MREstimate:
    """Summary of an MR evidence estimate.

    The optional intercept fields are populated only for MR-Egger.  The F
    statistics summarize instrument relevance and are diagnostic rather than
    proof that the other MR assumptions hold.
    """

    method: str
    beta: float
    se: float
    p_value: float
    ci_lower: float
    ci_upper: float
    n_instruments: int
    mean_f_statistic: float
    min_f_statistic: float
    q_statistic: float | None = None
    q_p_value: float | None = None
    intercept: float | None = None
    intercept_se: float | None = None
    intercept_p_value: float | None = None

    @property
    def estimate(self) -> float:
        """Alias for :attr:`beta`."""

        return self.beta

    @property
    def standard_error(self) -> float:
        """Alias for :attr:`se`."""

        return self.se

    @property
    def f_statistic(self) -> float:
        """Single-instrument F, or the mean F for multiple instruments."""

        return self.mean_f_statistic


def harmonize_alleles(
    exposure_effect_allele: str,
    exposure_other_allele: str,
    outcome_effect_allele: str,
    outcome_other_allele: str,
    outcome_beta: float,
    exposure_eaf: float | None = None,
    outcome_eaf: float | None = None,
    *,
    eaf_tolerance: float = 0.10,
) -> HarmonizationResult:
    """Align an outcome beta to the exposure effect allele.

    Direct, swapped, strand-complemented, and strand-complemented/swapped
    allele pairs are handled.  Palindromic A/T and C/G variants are retained
    only when both effect-allele frequencies identify one orientation within
    ``eaf_tolerance``; this avoids silently assigning the wrong sign near a
    frequency of 0.5.

    Status values are ``aligned``, ``flipped``, ``strand_aligned``,
    ``strand_flipped``, ``palindromic_aligned``,
    ``palindromic_flipped``, ``ambiguous_palindromic``, or ``incompatible``.
    """

    ea_x, oa_x = _validate_allele_pair(
        exposure_effect_allele, exposure_other_allele, "exposure"
    )
    ea_y, oa_y = _validate_allele_pair(
        outcome_effect_allele, outcome_other_allele, "outcome"
    )
    beta_y = _finite_scalar(outcome_beta, "outcome_beta")
    tolerance = _finite_scalar(eaf_tolerance, "eaf_tolerance")
    if not 0.0 < tolerance < 0.5:
        raise ValueError("eaf_tolerance must be strictly between 0 and 0.5")

    exposure_pair = (ea_x, oa_x)
    outcome_pair = (ea_y, oa_y)
    reversed_exposure = (oa_x, ea_x)
    complemented_outcome = (_COMPLEMENT[ea_y], _COMPLEMENT[oa_y])

    direct = outcome_pair == exposure_pair
    swapped = outcome_pair == reversed_exposure
    strand_direct = complemented_outcome == exposure_pair
    strand_swapped = complemented_outcome == reversed_exposure

    if not any((direct, swapped, strand_direct, strand_swapped)):
        return HarmonizationResult(float("nan"), "incompatible")

    is_palindromic = frozenset(exposure_pair) in _PALINDROMIC
    if is_palindromic:
        if exposure_eaf is None or outcome_eaf is None:
            return HarmonizationResult(float("nan"), "ambiguous_palindromic")
        eaf_x = _validate_eaf(exposure_eaf, "exposure_eaf")
        eaf_y = _validate_eaf(outcome_eaf, "outcome_eaf")

        no_flip_distance = abs(eaf_x - eaf_y)
        flip_distance = abs(eaf_x - (1.0 - eaf_y))
        no_flip_supported = no_flip_distance <= tolerance
        flip_supported = flip_distance <= tolerance

        # Near 0.5, both orientations can satisfy the frequency tolerance.
        # Likewise, frequencies supporting neither orientation are unsafe.
        if no_flip_supported == flip_supported:
            return HarmonizationResult(float("nan"), "ambiguous_palindromic")
        if no_flip_supported:
            return HarmonizationResult(beta_y, "palindromic_aligned")
        return HarmonizationResult(-beta_y, "palindromic_flipped")

    if direct:
        return HarmonizationResult(beta_y, "aligned")
    if swapped:
        return HarmonizationResult(-beta_y, "flipped")
    if strand_direct:
        return HarmonizationResult(beta_y, "strand_aligned")
    return HarmonizationResult(-beta_y, "strand_flipped")


def harmonize_dataframe(
    frame: pd.DataFrame,
    *,
    exposure_effect_allele_col: str = "exposure_effect_allele",
    exposure_other_allele_col: str = "exposure_other_allele",
    outcome_effect_allele_col: str = "outcome_effect_allele",
    outcome_other_allele_col: str = "outcome_other_allele",
    outcome_beta_col: str = "beta_y",
    exposure_eaf_col: str | None = "exposure_eaf",
    outcome_eaf_col: str | None = "outcome_eaf",
    eaf_tolerance: float = 0.10,
) -> pd.DataFrame:
    """Return a copy with signed ``beta_y_aligned`` and harmonization status.

    Rows are not silently removed.  Consumers can filter the returned
    ``harmonization_retained`` flag while preserving an audit trail of dropped
    palindromic or incompatible records.
    """

    if not isinstance(frame, pd.DataFrame):
        raise TypeError("frame must be a pandas DataFrame")
    required = {
        exposure_effect_allele_col,
        exposure_other_allele_col,
        outcome_effect_allele_col,
        outcome_other_allele_col,
        outcome_beta_col,
    }
    if exposure_eaf_col is not None:
        required.add(exposure_eaf_col)
    if outcome_eaf_col is not None:
        required.add(outcome_eaf_col)
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"frame is missing required columns: {missing}")

    results: list[HarmonizationResult] = []
    for row_index in range(len(frame)):
        row = frame.iloc[row_index]
        exposure_eaf = (
            None
            if exposure_eaf_col is None or pd.isna(row[exposure_eaf_col])
            else row[exposure_eaf_col]
        )
        outcome_eaf = (
            None
            if outcome_eaf_col is None or pd.isna(row[outcome_eaf_col])
            else row[outcome_eaf_col]
        )
        results.append(
            harmonize_alleles(
                row[exposure_effect_allele_col],
                row[exposure_other_allele_col],
                row[outcome_effect_allele_col],
                row[outcome_other_allele_col],
                row[outcome_beta_col],
                exposure_eaf,
                outcome_eaf,
                eaf_tolerance=eaf_tolerance,
            )
        )

    output = frame.copy()
    output["beta_y_aligned"] = [result.aligned_beta for result in results]
    output["harmonization_status"] = [result.status for result in results]
    output["harmonization_retained"] = [result.retained for result in results]
    return output


def wald_ratio(
    beta_x: float,
    se_x: float,
    beta_y: float,
    se_y: float,
    *,
    covariance: float = 0.0,
    confidence_level: float = 0.95,
) -> MREstimate:
    """Estimate a single-instrument Wald ratio using the full delta SE.

    The variance includes uncertainty in both ``beta_x`` and ``beta_y`` and an
    optional covariance term::

        Var(Y/X) = Var(Y)/X^2 + Y^2 Var(X)/X^4
                   - 2 Y Cov(X,Y)/X^3.

    Independent two-sample estimates normally use ``covariance=0``.
    """

    bx = _finite_scalar(beta_x, "beta_x")
    by = _finite_scalar(beta_y, "beta_y")
    sx = _positive_scalar(se_x, "se_x")
    sy = _positive_scalar(se_y, "se_y")
    cov = _finite_scalar(covariance, "covariance")
    _validate_covariance(cov, sx, sy)
    if bx == 0.0:
        raise ValueError("beta_x must be non-zero for a Wald ratio")

    # Express the delta variance as g' Sigma g using standardized derivative
    # terms.  This avoids raw x**4 calculations, which can underflow for valid
    # but extremely small floating-point inputs.
    with np.errstate(over="ignore", under="ignore", divide="ignore", invalid="ignore"):
        beta_value = np.divide(np.float64(by), np.float64(bx))
        outcome_term = np.divide(np.float64(sy), np.float64(bx))
        exposure_term = (
            -np.divide(np.float64(by), np.float64(bx))
            * np.divide(np.float64(sx), np.float64(bx))
        )
        correlation = (
            np.float64(0.0)
            if cov == 0.0
            else (
                np.divide(
                    np.divide(np.float64(cov), np.float64(sx)),
                    np.float64(sy),
                )
            )
        )
        variance_value = (
            np.square(outcome_term)
            + np.square(exposure_term)
            + 2.0 * outcome_term * exposure_term * correlation
        )
    if not np.isfinite(beta_value):
        raise ValueError("Wald ratio is non-finite; beta_x is too close to zero")
    if not np.isfinite(variance_value):
        raise ValueError(
            "delta-method variance is non-finite; instrument scale is unstable"
        )
    beta = float(beta_value)
    variance = float(variance_value)
    # A valid covariance matrix implies a non-negative delta variance.  Allow
    # only round-off-sized negative values before clipping.
    tolerance = np.finfo(float).eps * max(
        1.0, float(np.square(outcome_term)), float(np.square(exposure_term))
    ) * 32.0
    if variance < -tolerance:
        raise ValueError("delta-method variance is negative; check covariance")
    se = float(np.sqrt(max(0.0, variance)))
    if se == 0.0:
        raise ValueError("delta-method standard error is zero")

    ci_lower, ci_upper, p_value = _normal_summary(
        beta, se, confidence_level
    )
    with np.errstate(over="ignore", under="ignore"):
        f_statistic = float(np.square(np.divide(np.float64(bx), np.float64(sx))))
    if not np.isfinite(f_statistic):
        raise ValueError("F statistic is non-finite; input scales are unstable")
    return MREstimate(
        method="Wald ratio (full delta method)",
        beta=float(beta),
        se=se,
        p_value=p_value,
        ci_lower=ci_lower,
        ci_upper=ci_upper,
        n_instruments=1,
        mean_f_statistic=float(f_statistic),
        min_f_statistic=float(f_statistic),
    )


def inverse_variance_weighted(
    beta_x: ArrayLike,
    se_x: ArrayLike,
    beta_y: ArrayLike,
    se_y: ArrayLike,
    *,
    covariance: ArrayLike | float = 0.0,
    random_effects: bool = False,
    confidence_level: float = 0.95,
) -> MREstimate:
    """Combine at least two independent ratio estimates by inverse variance.

    Ratio variances use the same full delta method as :func:`wald_ratio`.
    With ``random_effects=True``, a multiplicative over-dispersion factor
    ``max(1, Q/(K-1))`` inflates the fixed-effect standard error.
    """

    bx, sx, by, sy, cov = _mr_arrays(beta_x, se_x, beta_y, se_y, covariance)
    if bx.size < 2:
        raise ValueError("IVW requires at least two independent instruments")
    ratios, variances = _ratio_estimates(bx, sx, by, sy, cov)
    weights = 1.0 / variances
    beta = float(np.dot(weights, ratios) / np.sum(weights))
    q_statistic = float(np.dot(weights, (ratios - beta) ** 2))
    q_p_value = float(chi2.sf(q_statistic, df=bx.size - 1))
    dispersion = (
        max(1.0, q_statistic / (bx.size - 1)) if random_effects else 1.0
    )
    se = float(np.sqrt(dispersion / np.sum(weights)))
    ci_lower, ci_upper, p_value = _normal_summary(
        beta, se, confidence_level
    )
    f_values = (bx / sx) ** 2
    return MREstimate(
        method=(
            "Inverse-variance weighted (multiplicative random effects)"
            if random_effects
            else "Inverse-variance weighted (fixed effect)"
        ),
        beta=beta,
        se=se,
        p_value=p_value,
        ci_lower=ci_lower,
        ci_upper=ci_upper,
        n_instruments=int(bx.size),
        mean_f_statistic=float(np.mean(f_values)),
        min_f_statistic=float(np.min(f_values)),
        q_statistic=q_statistic,
        q_p_value=q_p_value,
    )


def weighted_median(
    beta_x: ArrayLike,
    se_x: ArrayLike,
    beta_y: ArrayLike,
    se_y: ArrayLike,
    *,
    covariance: ArrayLike | float = 0.0,
    n_bootstrap: int = 5_000,
    random_seed: int = 0,
    confidence_level: float = 0.95,
) -> MREstimate:
    """Estimate the inverse-variance weighted median ratio.

    At least three independent instruments are required.  The standard error
    is obtained by a deterministic parametric bootstrap that propagates the
    supplied exposure and outcome uncertainty (and covariance).  The weighted
    median is robust when more than half of the total weight comes from valid
    instruments, conditional on the usual MR assumptions.
    """

    bx, sx, by, sy, cov = _mr_arrays(beta_x, se_x, beta_y, se_y, covariance)
    if bx.size < 3:
        raise ValueError(
            "weighted median requires at least three independent instruments"
        )
    if isinstance(n_bootstrap, (bool, np.bool_)) or not isinstance(
        n_bootstrap, (int, np.integer)
    ):
        raise TypeError("n_bootstrap must be an integer")
    if n_bootstrap < 2:
        raise ValueError("n_bootstrap must be at least 2")
    if isinstance(random_seed, (bool, np.bool_)) or not isinstance(
        random_seed, (int, np.integer)
    ):
        raise TypeError("random_seed must be an integer")
    if random_seed < 0:
        raise ValueError("random_seed must be non-negative")

    ratios, variances = _ratio_estimates(bx, sx, by, sy, cov)
    weights = 1.0 / variances
    beta = _weighted_quantile_median(ratios, weights)

    rng = np.random.default_rng(int(random_seed))
    z_x = rng.standard_normal((int(n_bootstrap), bx.size))
    z_y = rng.standard_normal((int(n_bootstrap), bx.size))
    rho = cov / (sx * sy)
    simulated_x = bx + sx * z_x
    simulated_y = by + sy * (
        rho * z_x + np.sqrt(np.maximum(0.0, 1.0 - rho**2)) * z_y
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        simulated_ratios = simulated_y / simulated_x
    bootstrap_estimates = _row_weighted_medians(simulated_ratios, weights)
    finite_bootstrap = bootstrap_estimates[np.isfinite(bootstrap_estimates)]
    if finite_bootstrap.size < max(2, int(0.99 * n_bootstrap)):
        raise RuntimeError(
            "too many non-finite bootstrap ratios; instruments may be too weak"
        )
    se = float(np.std(finite_bootstrap, ddof=1))
    if not np.isfinite(se) or se <= 0.0:
        raise RuntimeError("weighted-median bootstrap standard error is invalid")
    ci_lower, ci_upper, p_value = _normal_summary(
        beta, se, confidence_level
    )
    f_values = (bx / sx) ** 2
    return MREstimate(
        method="Weighted median (parametric-bootstrap SE)",
        beta=float(beta),
        se=se,
        p_value=p_value,
        ci_lower=ci_lower,
        ci_upper=ci_upper,
        n_instruments=int(bx.size),
        mean_f_statistic=float(np.mean(f_values)),
        min_f_statistic=float(np.min(f_values)),
    )


def mr_egger(
    beta_x: ArrayLike,
    se_x: ArrayLike,
    beta_y: ArrayLike,
    se_y: ArrayLike,
    *,
    random_effects: bool = True,
    confidence_level: float = 0.95,
) -> MREstimate:
    """Fit a weighted MR-Egger regression with an unconstrained intercept.

    At least three independent instruments with non-identical exposure effects
    are required.  Each instrument is internally re-oriented to the
    exposure-increasing allele (both its exposure and outcome beta are flipped
    together when needed); this orientation is essential for an interpretable
    Egger intercept.  Outcome inverse-variance weights are used, so this is the
    conventional NOME estimator: material exposure-measurement error can cause
    regression dilution.  The intercept is a directional-pleiotropy diagnostic
    under the InSIDE assumption, not a general test that pleiotropy is absent.
    """

    bx, sx, by, sy, _ = _mr_arrays(beta_x, se_x, beta_y, se_y, 0.0)
    if bx.size < 3:
        raise ValueError("MR-Egger requires at least three independent instruments")

    orientation = np.sign(bx)
    bx = bx * orientation
    by = by * orientation

    design = np.column_stack((np.ones(bx.size, dtype=float), bx))
    weights = 1.0 / sy**2
    information = design.T @ (weights[:, None] * design)
    if np.linalg.matrix_rank(information) < 2:
        raise ValueError(
            "MR-Egger requires non-identical exposure effect estimates"
        )
    covariance_base = np.linalg.inv(information)
    coefficients = covariance_base @ (design.T @ (weights * by))
    residual = by - design @ coefficients
    q_statistic = float(np.dot(weights, residual**2))
    degrees_freedom = bx.size - 2
    q_p_value = float(chi2.sf(q_statistic, df=degrees_freedom))
    dispersion = (
        max(1.0, q_statistic / degrees_freedom) if random_effects else 1.0
    )
    coefficient_covariance = covariance_base * dispersion
    standard_errors = np.sqrt(np.diag(coefficient_covariance))

    intercept = float(coefficients[0])
    beta = float(coefficients[1])
    intercept_se = float(standard_errors[0])
    se = float(standard_errors[1])
    intercept_p = float(2.0 * norm.sf(abs(intercept / intercept_se)))
    ci_lower, ci_upper, p_value = _normal_summary(
        beta, se, confidence_level
    )
    f_values = (bx / sx) ** 2
    return MREstimate(
        method=(
            "MR-Egger (multiplicative random effects)"
            if random_effects
            else "MR-Egger (fixed effect)"
        ),
        beta=beta,
        se=se,
        p_value=p_value,
        ci_lower=ci_lower,
        ci_upper=ci_upper,
        n_instruments=int(bx.size),
        mean_f_statistic=float(np.mean(f_values)),
        min_f_statistic=float(np.min(f_values)),
        q_statistic=q_statistic,
        q_p_value=q_p_value,
        intercept=intercept,
        intercept_se=intercept_se,
        intercept_p_value=intercept_p,
    )


def bh_fdr(p_values: ArrayLike) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values, preserving missing positions.

    NaNs are excluded from the number of tested hypotheses and remain NaN in
    the output.  Infinite values and finite values outside ``[0, 1]`` are
    rejected.  A stable sort makes tied inputs deterministic.
    """

    try:
        values = np.asarray(p_values, dtype=float)
    except (TypeError, ValueError) as exc:
        raise TypeError("p_values must be numeric") from exc
    if values.ndim != 1:
        raise ValueError("p_values must be one-dimensional")
    if np.any(np.isinf(values)):
        raise ValueError("p_values must not contain infinities")
    finite = np.isfinite(values)
    if np.any((values[finite] < 0.0) | (values[finite] > 1.0)):
        raise ValueError("finite p_values must lie in [0, 1]")

    adjusted = np.full(values.shape, np.nan, dtype=float)
    observed = values[finite]
    if observed.size == 0:
        return adjusted
    order = np.argsort(observed, kind="mergesort")
    ranked = observed[order]
    raw_adjusted = ranked * observed.size / np.arange(1, observed.size + 1)
    monotone = np.minimum.accumulate(raw_adjusted[::-1])[::-1]
    restored = np.empty_like(monotone)
    restored[order] = np.clip(monotone, 0.0, 1.0)
    adjusted[finite] = restored
    return adjusted


# Concise, commonly used alias.
ivw = inverse_variance_weighted


def _validate_allele_pair(
    effect_allele: str, other_allele: str, label: str
) -> tuple[str, str]:
    if not isinstance(effect_allele, str) or not isinstance(other_allele, str):
        raise TypeError(f"{label} alleles must be strings")
    pair = (effect_allele.strip().upper(), other_allele.strip().upper())
    if pair[0] not in _COMPLEMENT or pair[1] not in _COMPLEMENT:
        raise ValueError(f"{label} alleles must each be one of A, C, G, T")
    if pair[0] == pair[1]:
        raise ValueError(f"{label} effect and other alleles must differ")
    return pair


def _validate_eaf(value: float, name: str) -> float:
    eaf = _finite_scalar(value, name)
    if not 0.0 <= eaf <= 1.0:
        raise ValueError(f"{name} must lie in [0, 1]")
    return eaf


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


def _validate_covariance(covariance: float, se_x: float, se_y: float) -> None:
    if covariance == 0.0:
        return
    with np.errstate(over="ignore", under="ignore", divide="ignore", invalid="ignore"):
        correlation = np.divide(
            np.divide(np.float64(covariance), np.float64(se_x)),
            np.float64(se_y),
        )
    if not np.isfinite(correlation) or abs(correlation) > (
        1.0 + 32.0 * np.finfo(float).eps
    ):
        raise ValueError("covariance is incompatible with se_x and se_y")


def _mr_arrays(
    beta_x: ArrayLike,
    se_x: ArrayLike,
    beta_y: ArrayLike,
    se_y: ArrayLike,
    covariance: ArrayLike | float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    arrays: list[np.ndarray] = []
    for values, name in (
        (beta_x, "beta_x"),
        (se_x, "se_x"),
        (beta_y, "beta_y"),
        (se_y, "se_y"),
    ):
        try:
            array = np.asarray(values, dtype=float)
        except (TypeError, ValueError) as exc:
            raise TypeError(f"{name} must be numeric") from exc
        if array.ndim != 1:
            raise ValueError(f"{name} must be one-dimensional")
        arrays.append(array)
    bx, sx, by, sy = arrays
    if bx.size == 0:
        raise ValueError("at least one instrument is required")
    if not (bx.shape == sx.shape == by.shape == sy.shape):
        raise ValueError("beta_x, se_x, beta_y, and se_y must have equal length")
    for array, name in ((bx, "beta_x"), (sx, "se_x"), (by, "beta_y"), (sy, "se_y")):
        if not np.all(np.isfinite(array)):
            raise ValueError(f"{name} must contain only finite values")
    if np.any(bx == 0.0):
        raise ValueError("beta_x values must be non-zero")
    if np.any(sx <= 0.0) or np.any(sy <= 0.0):
        raise ValueError("se_x and se_y values must be strictly positive")

    try:
        cov = np.asarray(covariance, dtype=float)
    except (TypeError, ValueError) as exc:
        raise TypeError("covariance must be numeric") from exc
    if cov.ndim == 0:
        cov = np.full(bx.shape, float(cov), dtype=float)
    elif cov.ndim != 1 or cov.shape != bx.shape:
        raise ValueError("covariance must be scalar or match the instrument arrays")
    if not np.all(np.isfinite(cov)):
        raise ValueError("covariance must contain only finite values")
    with np.errstate(over="ignore", under="ignore", divide="ignore", invalid="ignore"):
        correlations = np.divide(np.divide(cov, sx), sy)
    if np.any(~np.isfinite(correlations)) or np.any(
        np.abs(correlations) > (1.0 + 32.0 * np.finfo(float).eps)
    ):
        raise ValueError("covariance is incompatible with se_x and se_y")
    return bx, sx, by, sy, cov


def _ratio_estimates(
    beta_x: np.ndarray,
    se_x: np.ndarray,
    beta_y: np.ndarray,
    se_y: np.ndarray,
    covariance: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    with np.errstate(over="ignore", under="ignore", divide="ignore", invalid="ignore"):
        ratios = beta_y / beta_x
        outcome_terms = se_y / beta_x
        exposure_terms = -(beta_y / beta_x) * (se_x / beta_x)
        correlations = np.divide(np.divide(covariance, se_x), se_y)
        variances = (
            outcome_terms**2
            + exposure_terms**2
            + 2.0 * outcome_terms * exposure_terms * correlations
        )
    scale = np.maximum.reduce(
        (
            np.ones_like(beta_x),
            outcome_terms**2,
            exposure_terms**2,
        )
    )
    tolerance = 32.0 * np.finfo(float).eps * scale
    if np.any(variances < -tolerance):
        raise ValueError("delta-method variance is negative; check covariance")
    variances = np.maximum(variances, 0.0)
    if np.any(~np.isfinite(variances)) or np.any(variances <= 0.0):
        raise ValueError("ratio variances must be finite and strictly positive")
    return ratios, variances


def _validate_confidence_level(confidence_level: float) -> float:
    level = _finite_scalar(confidence_level, "confidence_level")
    if not 0.0 < level < 1.0:
        raise ValueError("confidence_level must lie strictly between 0 and 1")
    return level


def _normal_summary(
    estimate: float, standard_error: float, confidence_level: float
) -> tuple[float, float, float]:
    level = _validate_confidence_level(confidence_level)
    quantile = float(norm.ppf(0.5 + level / 2.0))
    lower = float(estimate - quantile * standard_error)
    upper = float(estimate + quantile * standard_error)
    p_value = float(2.0 * norm.sf(abs(estimate / standard_error)))
    return lower, upper, p_value


def _weighted_quantile_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values, kind="mergesort")
    ordered_values = values[order]
    ordered_weights = weights[order]
    cumulative = np.cumsum(ordered_weights)
    half = 0.5 * cumulative[-1]
    index = int(np.searchsorted(cumulative, half, side="left"))
    if cumulative[index] == half and index + 1 < ordered_values.size:
        return float((ordered_values[index] + ordered_values[index + 1]) / 2.0)
    return float(ordered_values[index])


def _row_weighted_medians(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    order = np.argsort(values, axis=1, kind="stable")
    ordered_values = np.take_along_axis(values, order, axis=1)
    ordered_weights = weights[order]
    cumulative = np.cumsum(ordered_weights, axis=1)
    half = 0.5 * float(np.sum(weights))
    indices = np.argmax(cumulative >= half, axis=1)
    medians = ordered_values[np.arange(values.shape[0]), indices]

    exact = cumulative[np.arange(values.shape[0]), indices] == half
    can_average = exact & (indices + 1 < values.shape[1])
    rows = np.flatnonzero(can_average)
    if rows.size:
        medians[rows] = (
            medians[rows] + ordered_values[rows, indices[rows] + 1]
        ) / 2.0
    return medians


__all__ = [
    "HarmonizationResult",
    "MREstimate",
    "bh_fdr",
    "harmonize_alleles",
    "harmonize_dataframe",
    "inverse_variance_weighted",
    "ivw",
    "mr_egger",
    "wald_ratio",
    "weighted_median",
]
