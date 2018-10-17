'''
synbiochem (c) University of Manchester 2018

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from itertools import combinations
import sys

from ete3 import Tree
import pandas as pd


def prune(tree, num_nodes):
    '''Prune tree to representative number of nodes.'''
    dists = pd.DataFrame(index=tree.get_leaves(),
                         columns=tree.get_leaves())

    for pair in combinations(tree.get_leaves(), r=2):
        dists[pair[0]][pair[1]] = pair[0].get_distance(pair[1])

    while len(tree) > max(1, num_nodes):
        leaf1 = dists.min(axis=0).idxmin()
        leaf2 = dists.min(axis=1).idxmin()
        prune_leaf = leaf2 if leaf1.dist > leaf2.dist else leaf1
        dists.drop(prune_leaf, axis=0, inplace=True)
        dists.drop(prune_leaf, axis=1, inplace=True)

        _prune_tree(prune_leaf)


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
    tree = Tree(args[0], format=1)
    # t.resolve_polytomy()
    #  t.get_common_ancestor(n1, n2)
    # len(t)
    prune(tree, int(args[1]))

    print tree


if __name__ == '__main__':
    main(sys.argv[1:])
