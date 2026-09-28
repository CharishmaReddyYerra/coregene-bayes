from pathlib import Path
import json

import pandas as pd
import pytest

from coregene_bayes.bayesv2g_adapter import load_bayesv2g_candidates, top_seed_genes
from coregene_bayes.provenance import sha256_file


def test_bayesv2g_contract_with_minimal_fixture(tmp_path):
    frame = pd.DataFrame(
        [
            {
                "studyLocusId": "locus-1",
                "geneId": "ENSG1",
                "gene_symbol": "GENE1",
                "posterior_mean": 0.8,
                "posterior_lower": 0.6,
                "posterior_upper": 0.9,
                "conditional_posterior_mean": 0.8,
                "unresolved_probability": 0.1,
                "rank": 1,
                "probability_scope": "fixture",
                "evidence_coverage": 1.0,
            }
        ]
    )
    table = tmp_path / "prioritized_genes.csv"
    frame.to_csv(table, index=False)
    (tmp_path / "run_manifest.json").write_text(json.dumps({"fixture": True}), encoding="utf-8")
    (tmp_path / "SHA256SUMS").write_text(
        f"{sha256_file(table)}  prioritized_genes.csv\n", encoding="utf-8"
    )

    candidates, manifest = load_bayesv2g_candidates(tmp_path)

    assert top_seed_genes(candidates)["geneId"].tolist() == ["ENSG1"]
    assert manifest == {"fixture": True}


def test_real_bayesv2g_contract_is_verified():
    directory = Path(__file__).resolve().parents[2] / "bayesv2g" / "results" / "ldl_demo"
    if not directory.exists():
        pytest.skip("optional sibling BayesV2G real-artifact integration test")
    candidates, manifest = load_bayesv2g_candidates(directory)
    assert len(candidates) == 26
    assert len(top_seed_genes(candidates)) == 10
    assert isinstance(manifest, dict)
