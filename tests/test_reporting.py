from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd
import pytest

from coregene_bayes.reporting import plot_directed_network


def test_plot_directed_network_writes_nonempty_png(tmp_path) -> None:
    edges = pd.DataFrame(
        {
            "source_gene_symbol": ["APOB", "SORT1", "APOB"],
            "target_gene_symbol": ["SORT1", "LDLR", "LDLR"],
            "harmonization_retained": [True, True, False],
            "edge_probability": [0.91, 0.72, 0.99],
            "mr_beta": [0.15, -0.20, 0.40],
            "mr_p_value": [1e-8, 2e-5, 1e-10],
        }
    )
    path = tmp_path / "network.png"

    plot_directed_network(edges, path)

    assert path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
    assert path.stat().st_size > 1_000
    assert not plt.get_fignums()


def test_plot_directed_network_validates_limit(tmp_path) -> None:
    with pytest.raises(ValueError, match="positive"):
        plot_directed_network(pd.DataFrame(), tmp_path / "network.png", limit=0)
