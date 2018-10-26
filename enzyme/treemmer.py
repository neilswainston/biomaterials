'''
synbiochem (c) University of Manchester 2018

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
# pylint: disable=invalid-name
import sys

from ete3 import Tree
from synbiochem.utils import seq_utils

import pandas as pd


def get_homologues(id_seqs, num_nodes, common_ancestors=False, keep=None,
                   frmt=1):
    '''Get homologues.'''
    clustal_file = 'clustal.fasta'
    guidetree_file = 'tree.dnd'
    seq_utils.do_clustal(id_seqs, result_file=clustal_file,
                         guidetree_file=guidetree_file)

    return prune(guidetree_file, num_nodes, common_ancestors, keep, frmt)


def prune(filename, num_nodes, common_ancestors=False, keep=None, frmt=1):
    '''Prune tree to representative number of nodes.'''
    if keep is None:
        keep = []

    tree = Tree(filename, format=frmt)

    # Get subtree, based on common ancestry:
    sub_tree = _get_subtree(tree, common_ancestors, keep)

    # Form distance matrix:
    dists = _get_dists(sub_tree.get_leaves(), keep)

    # Iterate through distance matrix, deleting most redundant leaf:
    while len(sub_tree) > max(1, num_nodes):
        # Find most redundant leaf and delete:
        _remove_leaf(dists, keep)

    return sub_tree


def _get_subtree(tree, common_ancestors, keep):
    '''Get subtree.'''
    if common_ancestors and keep and len(keep) > 1:
        ancestors = list(keep)

        while len(ancestors) > 1:
            ancestor = tree.get_common_ancestor(ancestors[0], ancestors[1])
            ancestors = [ancestor] + ancestors[2:]

        return Tree(ancestors[0].write())

    return tree


def _get_dists(leaves, keep):
    '''Get distance matrix.'''
    dists = pd.DataFrame(index=leaves, columns=leaves)

    for leaf in leaves:
        _updated_dists(leaf, dists, keep)

    return dists


def _updated_dists(leaf, dists, keep):
    '''Update distance.'''
    for child in leaf.up.get_children():
        if child != leaf and child.is_leaf():
            # Order nodes by in keep, dist then name:
            nodes = sorted([leaf, child],
                           key=lambda node: (node.name in keep,
                                             node.dist,
                                             node.name))

            dists[nodes[0]][nodes[1]] = \
                float('inf') if all(node.name in keep for node in nodes) \
                else nodes[0].get_distance(nodes[1])


def _remove_leaf(dists, keep):
    '''Find most redundant leaf and delete.'''
    prune_leaf = dists.min(axis=0).idxmin()
    keep_leaf = dists.min(axis=1).idxmin()
    _prune_tree(prune_leaf)

    # Update distance matrix:
    _updated_dists(keep_leaf, dists, keep)
    dists.drop(prune_leaf, axis=0, inplace=True)
    dists.drop(prune_leaf, axis=1, inplace=True)


def _prune_tree(leaf):
    '''Prune leaf from tree.'''
    parent = leaf.up

    for child in parent.get_children():
        if child != leaf:
            child.dist += parent.dist

    leaf.detach()
    parent.delete()


def main(args):
    '''main method.'''
    id_seqs = seq_utils.read_fasta(args[0])
    print get_homologues(id_seqs, int(args[1]), args[2] != 'False', args[3:])


if __name__ == '__main__':
    main(sys.argv[1:])
