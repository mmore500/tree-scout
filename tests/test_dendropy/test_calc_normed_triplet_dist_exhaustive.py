import itertools as it
import math
import random
import typing

import alifedata_phyloinformatics_convert as apc
import networkx as nx
import pytest
import tqdist
from tqdm import tqdm

from tree_scout.dendropy._calc_normed_triplet_dist_exhaustive import (
    calc_normed_triplet_dist_exhaustive,
)


def all_possible_trees(n: int) -> typing.Iterator[nx.DiGraph]:
    for sequence in it.permutations(range(n), r=n - 2):
        undirected_graph = nx.from_prufer_sequence(sequence)
        for root_node in undirected_graph.nodes:
            res = nx.bfs_tree(undirected_graph, root_node).reverse()
            for node in res.nodes():
                res.nodes[node]["taxon_label"] = node
            yield res


def test_against_tqdist_exhaustive() -> None:
    for nx_tree in tqdm(all_possible_trees(6)):
        dendropy_tree = apc.RosettaTree(nx_tree).as_dendropy
        dendropy_tree.unassign_taxa(exclude_leaves=True)
        leaf_taxa = [leaf.taxon for leaf in dendropy_tree.leaf_node_iter()]
        if len(leaf_taxa) < 3:
            continue
        for leafing in it.permutations(leaf_taxa):
            comparison_tree = dendropy_tree.taxon_namespace_scoped_copy()
            for leaf, taxon in zip(comparison_tree.leaf_node_iter(), leafing):
                leaf.taxon = taxon
            assert math.isclose(
                tqdist.triplet_distance(
                    dendropy_tree.as_string(schema="newick").rstrip(),
                    comparison_tree.as_string(schema="newick").rstrip(),
                ),
                calc_normed_triplet_dist_exhaustive(
                    dendropy_tree,
                    comparison_tree,
                ),
            )


@pytest.mark.parametrize("tree_size", range(10, 30))
def test_identical_zero(tree_size: int) -> None:
    nx_tree = nx.random_tree(tree_size, create_using=nx.DiGraph).reverse()
    for node in nx_tree.nodes():
        nx_tree.nodes[node]["taxon_label"] = node
    dendropy_tree = apc.RosettaTree(nx_tree).as_dendropy

    leaf_taxa = [leaf.taxon for leaf in dendropy_tree.leaf_node_iter()]
    for taxon in leaf_taxa:
        assert taxon is not None
        assert taxon in dendropy_tree.taxon_namespace
    if len(leaf_taxa) < 3:
        return

    assert (
        calc_normed_triplet_dist_exhaustive(dendropy_tree, dendropy_tree) == 0
    )


@pytest.mark.parametrize("tree_size", range(10, 30))
def test_against_tqdist_sample(tree_size: int) -> None:
    nx_tree = nx.random_tree(
        tree_size,
        create_using=nx.DiGraph,
        seed=1,
    ).reverse()
    for node in nx_tree.nodes():
        nx_tree.nodes[node]["taxon_label"] = str(node)
    dendropy_tree = apc.RosettaTree(nx_tree).as_dendropy
    dendropy_tree.unassign_taxa(exclude_leaves=True)
    dendropy_tree.suppress_unifurcations()

    leaf_taxa = [leaf.taxon for leaf in dendropy_tree.leaf_node_iter()]
    if len(leaf_taxa) < 3:
        return
    comparison_tree = dendropy_tree.taxon_namespace_scoped_copy()
    comparison_tree.shuffle_taxa(rng=random.Random(1))

    # tqdist.triplet_distance must come first due to bug in dendropy
    # see https://github.com/jeetsukumaran/DendroPy/issues/152
    assert math.isclose(
        tqdist.triplet_distance(
            dendropy_tree.as_string(schema="newick").strip(),
            comparison_tree.as_string(schema="newick").strip(),
        ),
        calc_normed_triplet_dist_exhaustive(dendropy_tree, comparison_tree),
    )
