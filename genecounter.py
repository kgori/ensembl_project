#!/usr/bin/python
# -*- coding: utf-8 -*-

import ete2
import copy
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


class TreeAlignmentError(Exception):

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class NodeLabellingError(Exception):

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


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

    def identify(self):
        if not self.children:
            self.leaf = True
        else:
            if not self.up:
                self.root = True
            if hasattr(self, 'evoltype'):
                if self.evoltype == 'D':
                    self.duplication = True
                else:
                    self.speciation = True
            elif hasattr(self, 'D'):
                if self.D == 'Y':
                    self.duplication = True
                elif self.D == 'N':
                    self.speciation = True
            else:
                self.speciation = True

    def _just_gains(self, species_tree):
        self.identify()
        species_tree.identify()
        root = species_tree._find_mrca(self._taxon_set)
        root._justcount += 1
        for node in self.traverse():
            node.identify()
            if node.duplication:
                mrca = species_tree._find_mrca(node._taxon_set)
                mrca._justgains += 1

    def _get_species_subtree(self, species_tree, verbose=False):

        # define terms

        gtp = self  # gtp = GeneTree pointer, i.e. pointer to a node on the gene tree
        stp = species_tree  # stp = SpeciesTree pointer
        gtp.identify()
        stp.identify()
        if verbose: 
            print gtp
            if gtp.duplication:
                print 'gtp is duplication'
            if gtp.speciation:
                print 'gtp is speciation'
            if gtp.root:
                print 'gtp is root'
            if gtp.leaf:
                print 'gtp is leaf'
            print stp
            if stp.leaf:
                print 'stp is leaf'
            if stp.root:
                print 'stp is root'
        for taxon in gtp._taxon_set:
            if taxon not in stp._taxon_set:
                try:
                    raise TreeAlignmentError('Taxon {0} is not found in species tree'.format(taxon))
                except TreeAlignmentError, e:
                    print e.value


        if not stp.leaf:
            stp.L = species_tree.children[0]
            stp.R = species_tree.children[1]
            if verbose:
                print 'stp.L taxon set =', stp.L._taxon_set
                print 'stp.R taxon set =', stp.R._taxon_set

        if gtp.leaf:
            if stp.leaf:  # this is a match
                assert gtp._taxon_set == stp._taxon_set
                gtp._species_subtree = stp
                if gtp.root:
                    stp._count += 1
                return
            elif gtp._taxon_set <= stp.L._taxon_set:
                if not gtp.root:
                    stp.R._losses += 1
                gtp._get_species_subtree(stp.L)  # there have been deletions in the gene tree, lose stp.R and progress to stp.L
                return
            elif gtp._taxon_set <= stp.R._taxon_set:
                if not gtp.root:
                    stp.L._losses += 1
                gtp._get_species_subtree(stp.R, verbose)  # there have been deletions in the gene tree, lose stp.L and progress to stp.R
                return
            else:
                raise TreeAlignmentError('Gene tree leaf doesn\'t match species tree'
                        )
        else:
            gtp.L = self.children[0]
            gtp.R = self.children[1]
            if verbose:
                print 'gtp.L taxon set =', gtp.L._taxon_set
                print 'gtp.R taxon set =', gtp.R._taxon_set

        if stp.leaf:
            if verbose:
                print 'STP is a leaf'
                print 'GeneTreeNode =', gtp
                print 'SpeciesTreeNode =', stp
            try:
                assert gtp.duplication
            except AssertionError:
                try:
                    raise NodeLabellingError('Thought this node should be labelled as a duplication'
                            )
                except NodeLabellingError, e:
                    print 'NodeLabellingError:'
                    print gtp, stp, e
                    gtp.ambiguous = True
            gtp._species_subtree = stp
            stp._gains += 1
            if gtp.root:
                stp._count += 1
            gtp.L._get_species_subtree(stp, verbose)
            gtp.R._get_species_subtree(stp, verbose)
            return
        elif gtp.L._taxon_set <= stp.L._taxon_set \
            and gtp.L._taxon_set <= stp.R._taxon_set:
            if verbose:
                print 'GTP.L is a subset of both STP.L and STP.R'
                print 'GeneTreeNode =', gtp
                print 'SpeciesTreeNode =', stp
            raise TreeAlignmentError('''
GTP.L is a subset of both STP.L and STP.R
This implies that STP.L and STP.R intersect, which isn't allowed
on the species tree
'''
                    )
        elif gtp.L._taxon_set <= stp.L._taxon_set \
            and gtp.R._taxon_set <= stp.L._taxon_set:
            if verbose:
                print 'GTP.L is a subset of STP.L and GTP.R is a subset of STP.L - LOSS AT STP.R'
                print 'GeneTreeNode =', gtp
                print 'SpeciesTreeNode =', stp
            if not gtp.root:
                stp.R._losses += 1
            gtp._get_species_subtree(stp.L, verbose)
            return
        elif gtp.L._taxon_set <= stp.L._taxon_set \
            and gtp.R._taxon_set <= stp.R._taxon_set:
            if verbose:
                print 'GTP.L is a subset of STP.L and GTP.R is a subset of STP.R - MATCH'
                if gtp.duplication: print 'gtp is duplication - GAIN'
                print 'GeneTreeNode =', gtp
                print 'SpeciesTreeNode =', stp
            gtp._species_subtree = stp
            if gtp.root:
                stp._count += 1
            if gtp.duplication:
                raise TreeAlignmentError('shouldn\'t be a duplication here...'
                        )
            gtp.L._get_species_subtree(stp.L, verbose)
            gtp.R._get_species_subtree(stp.R, verbose)
            return
        elif gtp.L._taxon_set <= stp.R._taxon_set \
            and gtp.R._taxon_set <= stp.L._taxon_set:
            if verbose:
                print 'GTP.L is a subset of STP.R and GTP.R is a subset of STP.L - REVERSE MATCH'
                if gtp.duplication: print 'gtp is duplication - GAIN'
                print 'GeneTreeNode =', gtp
                print 'SpeciesTreeNode =', stp
            (gtp.R, gtp.L) = (gtp.L, gtp.R)
            gtp.swap_children()
            gtp._species_subtree = stp
            if gtp.root:
                stp._count += 1
            if gtp.duplication:
                raise TreeAlignmentError('shouldn\'t be a duplication here...'
                        )
            gtp.L._get_species_subtree(stp.L, verbose)
            gtp.R._get_species_subtree(stp.R, verbose)
            return
        elif gtp.L._taxon_set <= stp.R._taxon_set \
            and gtp.R._taxon_set <= stp.R._taxon_set:
            if verbose:
                print 'GTP.L is a subset of STP.R and GTP.R is a subset of STP.R - LOSS AT STP.L'
                print 'GeneTreeNode =', gtp
                print 'SpeciesTreeNode =', stp
            if not gtp.root:
                stp.L._losses += 1
            gtp._get_species_subtree(stp.R, verbose)
            return
        elif gtp.R._taxon_set <= stp.R._taxon_set \
            and gtp.R._taxon_set <= stp.L._taxon_set:
            if verbose:
                print 'GTP.R is a subset of STP.L and of STP.R'
                print 'GeneTreeNode =', gtp
                print 'SpeciesTreeNode =', stp
            raise TreeAlignmentError('''
GTP.R is a subset of STP.L and of STP.R
This implies that STP.L and STP.R intersect, which isn't allowed
on the species tree
'''
                    )
        elif stp.L._taxon_set & gtp._taxon_set and stp.R._taxon_set \
            & gtp._taxon_set:
            try:
                assert gtp.duplication
            except AssertionError:
                try:
                    raise NodeLabellingError('''
The gene tree topology doesn\'t match the species tree.
Maybe this node should be marked as a duplication?
Continuing as if this node were a duplication...
''')
                except NodeLabellingError, e:
                    print 'NodeLabellingError:'
                    print gtp, stp, e.value
                    compatibles = [gtp.L._taxon_set & stp.L._taxon_set, gtp.R._taxon_set & stp.R._taxon_set]
                    conflicts = [gtp.L._taxon_set & stp.R._taxon_set, gtp.R._taxon_set & stp.L._taxon_set]
                    print 'conflicts (l,r)', conflicts
                    print 'compatibles (l,r)', compatibles
                    gtp.ambiguous = True                

            stp._gains += 1
            if gtp.root:
                stp._count += 1
            gtp._species_subtree = stp
            gtp.L._get_species_subtree(stp, verbose)
            gtp.R._get_species_subtree(stp, verbose)
            return
        else:
            raise TreeAlignmentError('Case fell through')
    

    def __init__(self, newick=None, format=0):

        # ## From ete2.Tree class definition

        self._children = []
        self._up = None
        self._dist = 1.0
        self._support = 1.0
        self.features = set([])
        self._taxon_set = None  # added
        self._species_subtree = None  # added
        self.root = False  # added
        self.duplication = False  # added
        self.speciation = False  # added
        self.ambiguous = False  # added
        self.leaf = False  # added

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

    def display(self):
        for n in self.traverse():
            nstyle = ete2.NodeStyle()
            if n.duplication:
                nstyle['fgcolor'] = 'red'
                nstyle['size'] = 8
            elif n.speciation:
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

    def identify(self):
        if not self.children:
            self.leaf = True
        else:
            if not self.up:
                self.root = True

    def __init__(self, newick=None, format=0):

        # ## From ete2.Tree class definition

        self._children = []
        self._up = None
        self._dist = 1.0
        self._support = 1.0
        self._count = 0  # added
        self._gains = 0  # added
        self._justgains = 0
        self._justcount = 0
        self._losses = 0  # added
        self._change = 0  # added
        self._events = []  # added
        self.leaf = False  # added
        self.root = False  # added
        self._taxon_set = None  # added

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

    def display(self):
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
            if n.ambiguous:
                nstyle['fgcolor'] = 'cyan'
                nstyle['size'] = 6
            elif n.duplication:
                nstyle['fgcolor'] = 'red'
                nstyle['size'] = 6
            elif n.speciation:
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
        if len(sys.argv)>1: ts.title.add_face(ete2.TextFace(sys.argv[2], fsize=20),
                          column=0)
        genetree.show(tree_style=ts)

    def _find_mrca(self, taxon_set):
        """
        Returns the most recent common ancestor of the taxon set of the node
        """

        if len(taxon_set) > 1:
            return self.get_common_ancestor(taxon_set)
        else:
            return self.get_leaves_by_name(list(taxon_set)[0])[0]

    def _finalise_counts(self):
        self._count += self._gains - self._losses
        for node in self.iter_descendants('preorder'):
            node._count += node.up._count + node._gains - node._losses

