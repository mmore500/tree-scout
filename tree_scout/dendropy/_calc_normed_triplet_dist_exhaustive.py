import itertools as it
import statistics

import dendropy as dp
import numpy as np

from ._classify_triplet_topology import classify_triplet_topology


def calc_normed_triplet_dist_exhaustive(
    first: dp.Tree,
    second: dp.Tree,
) -> float:
    assert first.taxon_namespace is second.taxon_namespace
    taxa = [node.taxon for node in first.leaf_node_iter()]
    assert len(second.leaf_nodes()) == len(taxa)

    if len(taxa) < 3:
        return np.nan

    return 1.0 - statistics.mean(
        classify_triplet_topology(first, second, triplet)
        for triplet in it.combinations(taxa, r=3)
    )
