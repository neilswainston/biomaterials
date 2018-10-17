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
import random
import sys

from ete3 import Tree


def prune(tree, num_nodes):
    '''Prune tree to representative number of nodes.'''
    while len(tree) > max(3, num_nodes):
        neighbours = {}

        for leaf in tree.get_leaves():
            neighbours.update(_get_neighbours(leaf))

        prune_leaf = _get_prune_leaf(neighbours, tree)

        _prune_tree(prune_leaf, tree)

        # Purge the distance list of all pairs that have the pruned leaf:
        for key in [key for key in neighbours if prune_leaf in key]:
            del neighbours[key]


def _get_neighbours(leaf):
    '''Find neighbours of leaf.'''
    neighbours = {}
    parent = leaf.up
    flag = False
    sister_flag = False

    if parent.is_root():
        return neighbours
    else:
        # this for loop start from parent and climb up max two nodes,
        # if it finds leaves calculate the distances
        for n in range(0, len(parent.get_children())):
            # search at one node of distance
            if parent.children[n].is_leaf():
                if parent.children[n] != leaf:
                    dis = leaf.get_distance(parent.children[n])
                    neighbours[(leaf.name, parent.children[n].name)] = dis
                    flag = True
            elif not flag:
                temp_dlist = {}

                for nn in range(0, len(parent.children[n].get_children())):
                    if parent.children[n].children[nn].is_leaf():
                        dis = leaf.get_distance(
                            parent.children[n].children[nn])
                        temp_dlist[(leaf.name,
                                    parent.children[n].children[nn].name)] = dis
                        sister_flag = True

    # collect results at two nodes of distance only if there are no leaves
    # that are closer
    if sister_flag and not flag:
        neighbours.update(temp_dlist)

    if not flag:
        # this means that the leaf has no neighbors at one node of dist
        # therefore I climb the tree down towards the root of one more step and
        # look for leaves
        parent = parent.up
        multi_flag = False
        temp_dlist = {}
        # this for loop start from gran parent and climb up max one nodes, if
        # it finds leaves calculate the distances,
        for n in range(0, len(parent.get_children())):
            if parent.is_root():
                break
            if parent.children[n].is_leaf():
                dis = leaf.get_distance(parent.children[n])
                multi_flag = True
                temp_dlist[(leaf.name, parent.children[n].name)] = dis

        if multi_flag:					# this is to deal with polytomies
            neighbours.update(temp_dlist)

    return neighbours


def _get_prune_leaf(neighbours, t):
    '''parse the list with all neighbor pairs and distances, find the closest
    pair and select the leaf.'''
    min_val = min(neighbours.values())
    d_min = {}

    for k, v in neighbours.iteritems():
        if v == min_val:
            d_min.update({k: v})

    pair = random.choice(list(d_min))

    leaf1 = t.search_nodes(name=pair[0])[0]
    leaf2 = t.search_nodes(name=pair[1])[0]

    if leaf1.dist > leaf2.dist:
        leaf_to_prune = leaf2.name
        # leaf_to_keep = leaf1.name
    else:
        leaf_to_prune = leaf1.name
        # leaf_to_keep = leaf2.name

    # check if leaf is protected
    # if is_protected
        # leaf_to_prune = leaf_to_keep
        # if both leaves of the pair are protected => delete the pair from
        # neighbours and make another cycle
        # del neighbours[pair_unsplit]

    return leaf_to_prune


def _prune_tree(leaf_to_prune, tree):
    '''Prune leaf from tree.'''
    node = tree.search_nodes(name=leaf_to_prune)[0]
    parent = node.up

    if len(parent.get_children()) == 2:
        if parent.children[0] != node:
            parent.children[0].dist = parent.children[0].dist + parent.dist
        elif parent.children[1] != node:
            parent.children[1].dist = parent.children[1].dist + parent.dist

    node.detach()

    if len(parent.get_children()) == 1:
        # after pruning the remaining branch will be like this:
        # ---/---leaf_name. I delete useless node keeping the b length
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
