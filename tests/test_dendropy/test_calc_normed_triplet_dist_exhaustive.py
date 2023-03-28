import itertools as it
import math
import typing

import alifedata_phyloinformatics_convert as apc
import networkx as nx
import tqdist
from tqdm import tqdm

from tree_scout.dendropy._calc_normed_triplet_dist_exhaustive import (
    calc_normed_triplet_dist_exhaustive,
)


def all_possible_trees(n: int) -> typing.Iterator[nx.DiGraph]:
    for sequence in it.permutations(range(n), r=n - 2):
        undirected_graph = nx.from_prufer_sequence(sequence)
        for root_node in undirected_graph.nodes:
            yield nx.bfs_tree(undirected_graph, root_node).reverse()


def test_against_tqdist() -> None:
    for nx_tree in tqdm(all_possible_trees(6)):
        rosetta_tree = apc.RosettaTree(nx_tree)
        leaf_taxa = [
            leaf.taxon for leaf in rosetta_tree.as_dendropy.leaf_node_iter()
        ]
        if len(leaf_taxa) < 3:
            continue
        for leafing in it.permutations(leaf_taxa):
            comparison_tree = (
                rosetta_tree.as_dendropy.taxon_namespace_scoped_copy()
            )
            for leaf, taxon in zip(comparison_tree.leaf_node_iter(), leafing):
                leaf.taxon = taxon
                assert math.isclose(
                    tqdist.triplet_distance(
                        rosetta_tree.as_newick,
                        comparison_tree.as_string(schema="newick"),
                    ),
                    calc_normed_triplet_dist_exhaustive(
                        rosetta_tree.as_dendropy,
                        comparison_tree,
                    ),
                )
