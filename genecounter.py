#!/usr/bin/python
# -*- coding: utf-8 -*-

import ete2
import re
import sys
from ete2.parser.newick import read_newick, write_newick
from collections import defaultdict
try:
    from ete2.treeview.main import NodeStyle, _FaceAreas, \
        FaceContainer, FACE_POSITIONS
    from ete2.treeview.faces import Face
except ImportError:
    TREEVIEW = False
else:
    TREEVIEW = True


class GeneTree(ete2.PhyloTree):

    def parse_species_name(self):
        if not self.children:
            try:
                return self.name.split('_')[-1]
            except AttributeError:
                return self.name

    def get_species(self):
        for x in self.iter_leaves():
            yield x.parse_species_name()

    def _get_species_subtree(self, species_tree):
        """
        Recursive function for defining the subtree of the species tree associated
        with each node of the gene tree
        """

        # BASE CASE

        if not self.children:
            self._species_subtree = \
                species_tree.get_leaves_by_name(self.name)[0]
            return

        # VARIOUS OPTIONS FOR SETTING SPECIES SUBTREE

        while True:
            if self._species_subtree:  # subtree already set by being passed from parent speciation node
                break
            if self.type() == 'S':
                self._species_subtree = \
                    species_tree.get_common_ancestor(self._taxon_set)
                for child in self.children:
                    if child.type() == 'D':
                        for schild in self._species_subtree.children:
                            if child._taxon_set.issubset(schild._taxon_set):
                                child._species_subtree = schild
                break
            if not self.up:  # root node
                self._species_subtree = \
                    species_tree.get_common_ancestor(self._taxon_set)
                break
            if self.type() == 'D':
                self._species_subtree = self.up._species_subtree
                break

        # RECURSIVE CASE

        for child in self.children:
            child._get_species_subtree(species_tree)

    def __init__(self, newick=None, format=0):

        # ## From ete2.Tree class definition

        self._children = []
        self._up = None
        self._dist = 1.0
        self._support = 1.0
        self._taxon_set = None  # added
        self._species_subtree = None  # added
        self._index = None  # added

        self.features = set([])

        # Add basic features

        self.add_features(name='NoName')
        self.features.update(['dist', 'support'])

        # Initialize tree

        if newick is not None:
            read_newick(newick, root_node=self, format=format)

        if TREEVIEW:
            self._faces = _FaceAreas()
        for node in self.traverse():
            node._taxon_set = set(node.get_species())

        def nf(name):
            return name

        self.set_species_naming_function(nf)

    def type(self):
        try:
            if hasattr(self, 'evoltype'):
                return self.evoltype
            elif hasattr(self, 'D'):
                if self.children:
                    if self.D == 'Y':
                        return 'D'
                    elif self.D == 'N':
                        return 'S'
                else:
                    return 'L'
            else:
                return 'L'
        except AttributeError:
            return 'L'

    def display(self):
        for n in self.traverse():
            print n.type(),
            nstyle = ete2.NodeStyle()
            if n.type() == 'D':
                nstyle['fgcolor'] = 'red'
                nstyle['size'] = 8
            elif n.type() == 'S':
                nstyle['fgcolor'] = 'blue'
                nstyle['size'] = 8
            else:
                nstyle['fgcolor'] = 'lightblue'
            n.set_style(nstyle)
        self.show()


class SpeciesTree(ete2.Tree):

    def parse_species_name(self):
        if not self.children:
            try:
                return re.compile('(?<=_)[\w]+'
                                  ).search(self.name).group()
            except AttributeError:
                return self.name

    def get_species(self):
        for x in self.iter_leaves():
            yield x.parse_species_name()

    def __init__(self, newick=None, format=0):

        # ## From ete2.Tree class definition

        self._children = []
        self._up = None
        self._dist = 1.0
        self._support = 1.0
        self._count = 0  # added
        self._gains = 0  # added
        self._losses = 0  # added
        self._change = 0  # added
        self._events = []  # added
        self._taxon_set = None  # added
        self._already_deleted = []  # added

        self.features = set([])

        # Add basic features

        self.add_features(name='NoName')
        self.features.update(['dist', 'support'])

        # Initialize tree

        if newick is not None:
            read_newick(newick, root_node=self, format=format)

        if TREEVIEW:
            self._faces = _FaceAreas()
        for node in self.traverse():
            node._taxon_set = set(node.get_species())

    def traceback(self, node):
        while node.up and self != node:
            yield node
            node = node.up
        yield node

    def _find_mrca(self, taxon_set):
        """
        Returns the most recent common ancestor of the taxon set of the node
        """

        if len(taxon_set) > 1:
            return self.get_common_ancestor(taxon_set)
        else:
            return self.get_leaves_by_name(list(taxon_set)[0])[0]

    def display(self):
        for n in self.traverse():
            if hasattr(n, '_gains'):
                if n._gains > 0:
                    n.add_face(ete2.TextFace(n._gains, fgcolor='green'
                               ), column=0, position='branch-top')
            if hasattr(n, '_losses'):
                if n._losses:
                    n.add_face(ete2.TextFace(n._losses, fgcolor='red'),
                               column=0, position='branch-bottom')
        self.show()

    def display_with_gene_tree(self, genetree, scores=True):
        for n in self.traverse():
            if hasattr(n, '_gains'):
                if n._gains > 0:
                    n.add_face(ete2.TextFace('{0:+}'.format(n._gains),
                               fgcolor='green'), column=0,
                               position='branch-top')
            if hasattr(n, '_losses'):
                if n._losses:
                    n.add_face(ete2.TextFace('{0:+}'.format(-1
                               * n._losses), fgcolor='red'), column=0,
                               position='branch-bottom')
            if hasattr(n, '_count'):
                if n._count:
                    n.add_face(ete2.TextFace('{0}'.format(n._count),
                               fgcolor='blue'), column=0,
                               position='branch-right')
            nstyle = ete2.NodeStyle()
            if n._events:
                if 'Earliest evidence for this gene from init' \
                    in n._events:
                    nstyle['fgcolor'] = 'cyan'
            n.set_style(nstyle)
        for n in genetree.traverse():
            nstyle = ete2.NodeStyle()
            if n.type() == 'D':
                nstyle['fgcolor'] = 'red'
                nstyle['size'] = 6
            elif n.type() == 'S':
                nstyle['fgcolor'] = 'blue'
                nstyle['size'] = 6
            else:
                nstyle['fgcolor'] = 'lightblue'
            n.set_style(nstyle)

        def layout(node):
            species_ts = ete2.TreeStyle()
            species_ts.scale = 6
            species_ts.branch_vertical_margin = 1.5
            if not node.up:
                ete2.faces.add_face_to_node(ete2.TreeFace(self,
                        species_ts), node, 0, position='branch-top')

        ts = ete2.TreeStyle()
        ts.layout_fn = layout
        ts.title.add_face(ete2.TextFace(sys.argv[2], fsize=20),
                          column=0)
        genetree.show(tree_style=ts)

    def _find_losses(
        self,
        gene_tree_node,
        species_subtree=None,
        init=False,
        already_deleted=[],
        ):
        """
        inset = set of leaf nodes present in the genetree
        """

        inset = gene_tree_node._taxon_set
        if not species_subtree:
            species_subtree = self._find_mrca(inset)
        totalset = set(species_subtree.get_leaf_names())
        outset = totalset - inset
        traced = []
        untraced = defaultdict(list)
        deletions = []

        for inleaf in inset:
            inNode = self & inleaf
            for n in species_subtree.traceback(inNode):
                if n not in traced:
                    traced.append(n)
                else:
                    break

        for outleaf in outset:
            outNode = self & outleaf
            for n in species_subtree.traceback(outNode):
                if n not in traced:
                    untraced[outleaf].append(n)
                else:
                    break

        # each value in untraced is a list of all the untraced nodes
        # under each leaf in the outset. We only want to count the
        # earliest ones.

        for val in untraced.values():
            earliest = sorted(val, key=lambda x: self.get_distance(x,
                              topology_only=True))[0]
            if earliest not in deletions + already_deleted:
                deletions.append(earliest)

        if init:
            for node in deletions:
                node._events.append('Gene loss implied by init')
                node._losses += 1
        else:
            for node in deletions:
                node._events.append('Gene loss implied by duplication node'
                                    )
                node._losses += 1

        return deletions + already_deleted

    def _add_duplication_nodes(
        self,
        gene_tree_node,
        already_deleted=[],
        verbose=False,
        **kwargs
        ):

        # CASE 1 - root node

        if not gene_tree_node.up:
            if gene_tree_node.type() == 'D':
                if verbose:
                    print 'DUPLICATION NODE FOUND'
                if verbose:
                    print gene_tree_node, already_deleted
                mrca = self._find_mrca(gene_tree_node._taxon_set)
                mrca._gains += 1
                mrca._events.append('Gene duplication')
                for child in gene_tree_node.children:
                    self._add_duplication_nodes(child,
                            already_deleted=already_deleted)
            else:
                if verbose:
                    print 'ROOT NODE FOUND, RECURSING'
                for child in gene_tree_node.children:
                    self._add_duplication_nodes(child,
                            already_deleted=already_deleted)
        elif gene_tree_node.type() == 'D' and gene_tree_node.up.type() \
            == 'D':

        # CASE 2 - duplication following duplication                                    # This case checks for speciation followed by loss,
                                                                                        # which would result in a missing speciation node

            if verbose:  # that would otherwise be observed between the
                print 'DUPLICATION FOLLOWING DUPLICATION', \
                    already_deleted
            mrca = self._find_mrca(gene_tree_node._taxon_set)  # duplication nodes
            mrca._gains += 1
            mrca._events.append('Gene duplication')

            # subtree = gene_tree_node._species_subtree

            subtree = self._find_mrca(gene_tree_node.up._taxon_set)
            already_deleted = self._find_losses(gene_tree_node,
                    subtree, already_deleted=already_deleted)
            for d in already_deleted:
                if d not in self._already_deleted:
                    if verbose:
                        print 'adding', d
                    self._already_deleted.append(d)
            for child in gene_tree_node.children:
                self._add_duplication_nodes(child,
                        already_deleted=already_deleted)
        elif gene_tree_node.type() == 'D':

        # CASE 3 - duplication                                                          # Record gain at mrca, then recurse on child nodes

            if verbose:
                print 'DUPLICATION NODE FOUND'
            if verbose:
                print gene_tree_node, already_deleted
            mrca = self._find_mrca(gene_tree_node._taxon_set)
            mrca._gains += 1
            mrca._events.append('Gene duplication')
            for child in gene_tree_node.children:
                self._add_duplication_nodes(child,
                        already_deleted=already_deleted)
        elif gene_tree_node.type() != 'D' and gene_tree_node.up.type() \
            == 'D':

        # CASE 4 - speciation following duplication                                     # Using mrca found at parent duplication node,
                                                                                        # record deletions under this node in already_

            if verbose:  # -deleted list and _already_deleted global,
                print 'SPECIATION FOLLOWING DUPLICATION', \
                    already_deleted
            subtree = self._find_mrca(gene_tree_node.up._taxon_set)  # then recurse on child nodes
            already_deleted = self._find_losses(gene_tree_node,
                    subtree, already_deleted=already_deleted)
            for d in already_deleted:
                if d not in self._already_deleted:
                    if verbose:
                        print 'adding', d
                    self._already_deleted.append(d)
            if gene_tree_node.children:
                for child in gene_tree_node.children:
                    self._add_duplication_nodes(child,
                            already_deleted=already_deleted)
            else:
                return
        elif not gene_tree_node.children:

        # CASE 5 - leaf node                                                            # base case, end recursion

            if verbose:
                print 'LEAF NODE FOUND, RETURNING'
            if verbose:
                print gene_tree_node, already_deleted
            return
        else:

        # CASE 6 - anything else
                                                                                        # Just recurse on child nodes

            if verbose:
                print 'SPECIATION NODE FOUND', already_deleted
            if verbose:
                print gene_tree_node
            for child in gene_tree_node.children:
                self._add_duplication_nodes(child,
                        already_deleted=already_deleted)

    def _add_root_node(self, gene_tree_node):

        # if not gene_tree_node.type() == 'D':
        # print 'ADDING ROOT NODE'

        mrca = self._find_mrca(gene_tree_node._taxon_set)
        mrca._count += 1
        mrca._events.append('Earliest evidence for this gene from init')
        self._find_losses(gene_tree_node, init=True,
                          already_deleted=self._already_deleted)

    def _finalise_counts(self):
        self._count += self._gains - self._losses
        for node in self.iter_descendants('preorder'):
            node._count += node.up._count + node._gains - node._losses

    def add_gene_tree(self, gene_tree):

        # gene_tree._get_species_subtree(self)

        self._add_duplication_nodes(gene_tree)
        self._add_root_node(gene_tree)
        self._finalise_counts()
        self._already_deleted = []
