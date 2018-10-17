'''
# Treemmer

# Copyright 2018 Fabrizio Menardo


#   This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Dependencies:
# You have to install ete3 http://etetoolkit.org/
# and joblib https://pythonhosted.org/joblib/ to run Treemmer

# If you use Treemmer for your research, please cite:
# Treemmer: a tool to reduce large phylogenetic datasets with minimal loss
# of diversity. Menardo et. al (2018),BMC Bioinformatics 19:164.
# https://doi.org/10.1186/s12859-018-2164-8
'''
# pylint: disable=invalid-name
# pylint: disable=too-many-arguments
# pylint: disable=too-many-branches
# pylint: disable=too-many-locals
# pylint: disable=too-many-nested-blocks
# pylint: disable=too-many-statements
# from itertools import combinations
import operator
import sys

from ete3 import Tree


def prune(tree, num_nodes):
    '''Prune tree to representative number of nodes.'''
    while len(tree) > max(3, num_nodes):
        neighbours = []

        # for pair in combinations(tree.get_leaves(), r=2):
        #    neighbours.append([pair, pair[0].get_distance(pair[1])])

        for leaf in tree.get_leaves():
            neighbours.extend(_get_neighbours(leaf))

        prune_leaf = _get_prune_leaf(neighbours)
        _prune_tree(prune_leaf)


def _get_neighbours(leaf):
    '''Find neighbours of leaf.'''
    neighbours = []
    parent = leaf.up
    flag = False

    if parent.is_root():
        return neighbours
    else:
        # this for loop start from parent and climb up max two nodes,
        # if it finds leaves calculate the distances
        for child in parent.get_children():
            # search at one node of distance
            if child.is_leaf():
                if child != leaf:
                    dis = leaf.get_distance(child)
                    neighbours.append([leaf, child, dis])
                    flag = True
            elif not flag:
                for grandchild in child.get_children():
                    if grandchild.is_leaf():
                        dis = leaf.get_distance(grandchild)
                        neighbours.append([leaf, grandchild, dis])

    if not flag:
        # this means that the leaf has no neighbors at one node of dist
        # therefore I climb the tree down towards the root of one more step and
        # look for leaves
        parent = parent.up

        # this for loop start from gran parent and climb up max one nodes, if
        # it finds leaves calculate the distances,
        for n in range(0, len(parent.get_children())):
            if parent.is_root():
                break
            if parent.children[n].is_leaf():
                dis = leaf.get_distance(parent.children[n])
                neighbours.append([leaf, parent.children[n], dis])

    return neighbours


def _get_prune_leaf(neighbours):
    '''parse the list with all neighbor pairs and distances, find the closest
    pair and select the leaf.'''
    neighbours.sort(key=operator.itemgetter(2))

    pair = neighbours[0]

    return pair[1] if pair[0].dist > pair[1].dist else pair[0]


def _prune_tree(prune_leaf):
    '''Prune leaf from tree.'''
    parent = prune_leaf.up

    for child in parent.get_children():
        if child != prune_leaf:
            child.dist += parent.dist

    prune_leaf.detach()

    # after pruning the remaining branch will be like this:
    # ---/---leaf_name.
    # delete useless node keeping the b length
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
