import pandas as pd

from coregene_bayes.network import assign_physical_source_loci, network_coverage


def test_nearby_genes_share_physical_locus():
    genes = pd.DataFrame(
        {
            "gene_id": ["a", "b", "c", "d"],
            "gene_chromosome": ["1", "1", "1", "2"],
            "gene_start": [100, 900_000, 3_000_000, 50],
            "gene_end": [200, 900_100, 3_000_100, 100],
        }
    )
    loci = assign_physical_source_loci(genes, 1_000_000)
    assert loci.iloc[0] == loci.iloc[1]
    assert loci.iloc[0] != loci.iloc[2]
    assert loci.iloc[0] != loci.iloc[3]


def test_network_coverage_keeps_missing_as_untested():
    edges = pd.DataFrame(
        {
            "target_gene_id": ["x", "x", "y"],
            "source_locus_id": ["l1", "l2", "l1"],
            "harmonization_status": ["aligned", "aligned", "incompatible"],
            "harmonization_retained": [True, True, False],
        }
    )
    result = network_coverage(edges, {"l1", "l2", "l3", "l4"}, pd.Series(["x", "y", "z"]))
    observed = result.set_index("gene_id")
    assert observed.loc["x", "network_coverage"] == 0.5
    assert observed.loc["y", "tested_source_loci"] == 0
    assert observed.loc["z", "tested_source_loci"] == 0

