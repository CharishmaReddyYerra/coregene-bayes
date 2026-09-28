import numpy as np
import pytest

from coregene_bayes.coloc import coloc_abf, wakefield_log_abf


def test_coloc_prefers_shared_signal_when_peaks_match():
    beta_one = np.array([0.02, 0.75, 0.01])
    beta_two = np.array([0.01, -0.70, 0.02])
    se = np.repeat(0.05, 3)
    result = coloc_abf(beta_one, se, beta_two, se)
    assert result.pp_h4 > 0.95
    assert result.top_shared_index == 1
    assert sum(result.posterior.values()) == pytest.approx(1.0)


def test_coloc_prefers_distinct_signals_when_peaks_differ():
    se = np.repeat(0.05, 3)
    result = coloc_abf(np.array([0.8, 0.01, 0.01]), se, np.array([0.01, 0.8, 0.01]), se)
    assert result.posterior["H3"] > result.posterior["H4"]


def test_wakefield_rejects_bad_se():
    with pytest.raises(ValueError, match="standard errors"):
        wakefield_log_abf(np.array([1.0]), np.array([0.0]))

