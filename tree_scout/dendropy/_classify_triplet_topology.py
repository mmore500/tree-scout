import typing

import dendropy as dp


def classify_triplet_topology(
    tree: dp.Tree, triplet: typing.Sequence[dp.Taxon]
) -> typing.Optional[float]:
    first, second, third = triplet
    assert first in tree.taxon_namespace
    assert second in tree.taxon_namespace
    assert third in tree.taxon_namespace

    triplet_mrca = tree.mrca(taxa=triplet)
    mrca1 = tree.mrca(taxa=(first, second), start_node=triplet_mrca)
    mrca2 = tree.mrca(taxa=(first, third), start_node=triplet_mrca)
    mrca3 = tree.mrca(taxa=(second, third), start_node=triplet_mrca)
    if id(mrca1) == id(mrca2) == id(mrca3):
        return None
    elif id(mrca1) == id(mrca2):
        return 1
    elif id(mrca1) == id(mrca3):
        return 2
    elif id(mrca2) == id(mrca3):
        return 3
    else:
        assert False
