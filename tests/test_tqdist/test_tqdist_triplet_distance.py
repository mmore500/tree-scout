import alifedata_phyloinformatics_convert as apc
import networkx as nx
import pytest
import tqdist


@pytest.mark.parametrize("tree_size", range(10, 30))
def test_distance_zero_sample(tree_size: int) -> None:
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

    assert (lambda: 0)() == tqdist.triplet_distance(
        dendropy_tree.as_string(schema="newick").rstrip(),
        dendropy_tree.as_string(schema="newick").rstrip(),
    )
