import json

import pytest

from coregene_bayes.config import load_config


def test_default_config_loads():
    config = load_config()
    assert config["disease_study_id"] == "GCST90132314"
    assert config["core_thresholds"] == [2, 3, 5]


def test_invalid_probability_rejected(tmp_path):
    source = load_config()
    source.pop("config_path")
    source["edge_prior_inclusion"] = 1.2
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(source), encoding="utf-8")
    with pytest.raises(ValueError, match="edge_prior_inclusion"):
        load_config(path)
