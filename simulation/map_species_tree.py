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
            print 'GTP=', gtp
        if gtp.root:
            if verbose:
                print 'GTP is root'
        if gtp.speciation:
            if verbose:
                print 'GTP is speciation'
        if gtp.duplication:
            if verbose:
                print 'GTP is duplication'
        if gtp.leaf:
            if verbose:
                print 'GTP is leaf'
        if verbose:
            print 'STP=', stp
        if stp.leaf:
            if verbose:
                print 'STP is leaf'

        if not stp.leaf:
            if verbose:
                print 'stp is not a leaf'
            stp.L = species_tree.children[0]
            stp.R = species_tree.children[1]

        if gtp.leaf:
            if verbose:
                print 'gtp is a leaf'
            if stp.leaf:
                if verbose:
                    print 'stp is a leaf as well - MATCH'
                gtp._species_subtree = stp
                if gtp.root:
                    stp._count += 1
                    print 'l:124=root count +1', gtp,stp
                if verbose:
                    print 'gtp._species subtree set to stp'
                    print 'count, gains, losses=',stp._count, stp._gains, stp._losses
                if verbose:
                    print gtp._species_subtree
                return
            elif gtp._taxon_set <= stp.L._taxon_set:
                if verbose:
                    print 'gtp is a subset of stp.L, recursing...'
                if not gtp.root:
                    stp.R._losses += 1
                    print 'l:136 loss stp.R =count, gains, losses=',gtp,stp,stp.R._count, stp.R._gains, stp.R._losses
                gtp._get_species_subtree(stp.L, verbose)
                return
            elif gtp._taxon_set <= stp.R._taxon_set:
                if verbose:
                    print 'gtp is a subset of stp.R, recursing...'
                if not gtp.root:
                    stp.L._losses += 1
                    print 'l:144 loss stp.L =count, gains, losses=',gtp,stp,stp.L._count, stp.L._gains, stp.L._losses
                gtp._get_species_subtree(stp.R, verbose)
                return
            else:
                raise TreeAlignmentError('Gene tree leaf doesn\'t match species tree'
                        )
        else:
            if verbose:
                print 'gtp is not a leaf'
            gtp.L = self.children[0]
            gtp.R = self.children[1]

            if gtp.root:
                if verbose:
                    print 'gtp is root'
            if gtp.speciation:
                if verbose:
                    print 'gtp is a speciation'
            if gtp.duplication:
                if verbose:
                    print 'gtp is a duplication'

            if stp.leaf:
                if verbose:
                    print 'stp is leaf, but gtp isn\'t, implies duplication'
                try:
                    assert gtp.duplication
                except AssertionError:
                    try:
                        raise NodeLabellingError('Thought this node should be labelled as a duplication'
                                )
                    except NodeLabellingError, e:
                        print 'NodeLabellingError:'
                        print gtp, stp, e
                        pass
                if verbose:
                    print 'matching and recursing'
                gtp._species_subtree = stp
                stp._gains += 1
                print 'count, gains, losses=',stp._count, stp._gains, stp._losses
                if gtp.root:
                    stp._count += 1
                    print 'l:186, root count +1=count, gains, losses=',gtp,stp,stp._count, stp._gains, stp._losses
                gtp.L._get_species_subtree(stp, verbose)
                gtp.R._get_species_subtree(stp, verbose)
                return
            elif gtp.L._taxon_set <= stp.L._taxon_set \
                and gtp.L._taxon_set <= stp.R._taxon_set:
                if verbose:
                    print 'gtp.L is a subset of stp.L'
                if verbose:
                    print 'gtp.L is a subset of stp.R - ERROR'
                raise TreeAlignmentError('Possible Mislabelling of GeneTree nodes'
                        )
            elif gtp.L._taxon_set <= stp.L._taxon_set \
                and gtp.R._taxon_set <= stp.L._taxon_set:
                if verbose:
                    print 'gtp.L is a subset of stp.L'
                if verbose:
                    print 'gtp.R is a subset of stp.L'
                if verbose:
                    print 'advancing stp to stp.L, recursing....'
                if not gtp.root:
                    stp.R._losses += 1
                    print 'l:208 stp.R loss=count, gains, losses=',gtp,stp,stp.R._count, stp.R._gains, stp.R._losses
                gtp._get_species_subtree(stp.L, verbose)
                return
            elif gtp.L._taxon_set <= stp.L._taxon_set \
                and gtp.R._taxon_set <= stp.R._taxon_set:
                if verbose:
                    print 'gtp.L is a subset of stp.L'
                if verbose:
                    print 'gtp.R is a subset of stp.R'
                gtp._species_subtree = stp
                if verbose:
                    print 'MATCH. setting gtp._species_subtree to stp. recursing...'
                if gtp.root:
                    stp._count += 1
                    print 'l:222 - root gain =count, gains, losses=',gtp,stp,stp._count, stp._gains, stp._losses
                if gtp.duplication:
                    raise TreeAlignmentError('shouldn\'t be a duplication here...'
                            )
                gtp.L._get_species_subtree(stp.L, verbose)
                gtp.R._get_species_subtree(stp.R, verbose)
                return
            elif gtp.L._taxon_set <= stp.R._taxon_set \
                and gtp.R._taxon_set <= stp.L._taxon_set:
                if verbose:
                    print 'gtp.L is a subset of stp.R'
                if verbose:
                    print 'gtp.R is a subset of stp.L'
                (gtp.R, gtp.L) = (gtp.L, gtp.R)
                gtp.swap_children()
                gtp._species_subtree = stp
                if verbose:
                    print 'REVERSE MATCH. setting gtp._species_subtree to stp, recursing...'
                if gtp.root:
                    stp._count += 1
                    print 'l:242 - root gain=count, gains, losses=',gtp,stp,stp._count, stp._gains, stp._losses
                if gtp.duplication:
                    raise TreeAlignmentError('shouldn\'t be a duplication here...'
                            )
                gtp.L._get_species_subtree(stp.L, verbose)
                gtp.R._get_species_subtree(stp.R, verbose)
                return
            elif gtp.L._taxon_set <= stp.R._taxon_set \
                and gtp.R._taxon_set <= stp.R._taxon_set:
                if verbose:
                    print 'gtp.L is a subset of stp.R'
                if verbose:
                    print 'gtp.R is a subset of stp.R'
                if verbose:
                    print 'advancing pointer to stp.R, recursing...'
                if not gtp.root:
                    stp.L._losses += 1
                    print 'l:259, stp.L loss=count, gains, losses=',gtp,stp,stp.L._count, stp.L._gains, stp.L._losses
                gtp._get_species_subtree(stp.R, verbose)
                return
            elif gtp.R._taxon_set <= stp.R._taxon_set \
                and gtp.R._taxon_set <= stp.L._taxon_set:
                if verbose:
                    print 'gtp.R is a subset of stp.R and of stp.L - ERROR'
                raise TreeAlignmentError('Possible Mislabelling of GeneTree nodes'
                        )
            elif stp.L._taxon_set & gtp._taxon_set and stp.R._taxon_set \
                & gtp._taxon_set:

                if verbose:
                    print 'stp.L and stp.R intersect with union of gtp.L and gtp.R'
                if verbose:
                    print 'Implies duplication node, assertion check...'
                try:
                    assert gtp.duplication
                except AssertionError:
                    try:
                        raise NodeLabellingError('Thought this node should be labelled as a duplication'
                                )
                    except NodeLabellingError, e:
                        print 'NodeLabellingError:'
                        print gtp, stp, e
                        pass

                if verbose:
                    print 'assertion pass.'
                if verbose:
                    print 'MATCH, and gain'
                stp._gains += 1
                print 'l:291 gain=count, gains, losses=',gtp,stp,stp._count, stp._gains, stp._losses
                if gtp.root:
                    stp._count += 1
                    print 'l:294 root gain=count, gains, losses=',gtp,stp,stp._count, stp._gains, stp._losses
                gtp._species_subtree = stp
                gtp.L._get_species_subtree(stp, verbose)
                gtp.R._get_species_subtree(stp, verbose)
                return
            else:
                raise TreeAlignmentError('Case fell through')
        return

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
        self._justgains =0
        self._justcount =0
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
            if n.duplication:
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
        ts.title.add_face(ete2.TextFace(sys.argv[2], fsize=20),
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


# st = \
#     SpeciesTree('(SE001:60.971176,(SE002:54.097175,((SE003:8.863838,SE013:8.7271048):33.859374,(((SE004:40.303116,(SE007:35.957807,(SE008:27.55188,(SE010:23.561527,SE011:23.148886):2.2511886):0.10991615):3.753085):7.938699,SE006:35.196556):7.992514,((SE005:16.104527,(SE014:1.9441679,SE015:3.2796744):19.411262):15.863363,(SE009:19.907848,SE012:39.749202):5.4071946):19.687332):3.7807204):17.137093):7.2302568):0;'
#                 )
# gt = \
#     GeneTree('(SE001,(SE002,((((SE004,(SE007,(SE008,SE010))),SE006),((SE004,(SE007,(SE008,(SE008,SE010))[&&NHX:D=Y])),SE006))[&&NHX:D=Y],(((SE005,(SE014,SE015)),(SE009,SE012)),(SE005,(SE009,SE012)))[&&NHX:D=Y])));'
#              , format=1)
# gt._get_species_subtree(st)

# for node in gt.traverse():
#     if not node._species_subtree:
#         print 'NODE\n', node, node.type()
#         print 'SUBTREE\n', node._species_subtree

# st.display_with_gene_tree(gt)

st = SpeciesTree(open(sys.argv[1]).read(), format=1)
st2 = SpeciesTree(open(sys.argv[1]).read(), format=1)
gt = GeneTree(open(sys.argv[2]).read(), format=1)
if not '[&&NHX' in open(sys.argv[2]).read():
    gt.get_descendant_evol_events()

# st.display_with_gene_tree(gt)
# gt._just_gains(st)
gt._get_species_subtree(st, verbose=False)
st._finalise_counts()

# st.display_with_gene_tree(gt)
# for node in gt.traverse('preorder'):
#     new = node._species_subtree
#     old = node._species_subtree_OLD
#     mrca = st._find_mrca(node._taxon_set)
#     if new != old:
#         print 'NEW != OLD'
#         print 'NODE\n', node
#         print 'NEW SUBTREE\n', node._species_subtree
#         print 'OLD SUBTREE\n', node._species_subtree_OLD
#         print 'MRCA SUBTREE\n', st._find_mrca(node._taxon_set)
#         print
#         print
#         print
#     elif new != mrca:

#         # st.display_with_gene_tree(gt)

#         print 'NEW != MRCA'
#         print 'NODE\n', node
#         print 'NEW SUBTREE\n', node._species_subtree
#         print 'OLD SUBTREE\n', node._species_subtree_OLD
#         print 'MRCA SUBTREE\n', st._find_mrca(node._taxon_set)
#         print
#         print
#         print
#     elif mrca != old:

#         # st.display_with_gene_tree(gt)

#         print 'MRCA != OLD'
#         print 'NODE\n', node
#         print 'NEW SUBTREE\n', node._species_subtree
#         print 'OLD SUBTREE\n', node._species_subtree_OLD
#         print 'MRCA SUBTREE\n', st._find_mrca(node._taxon_set)
#         print
#         print
#         print

#         # st.display_with_gene_tree(gt)

for node in st.iter_leaves():
    print node.name, node._count, \
        list(gt.get_species()).count(node.name)
    try:
        assert node._count == list(gt.get_species()).count(node.name)
    except:
        print 'AssertionError'
        break

        # print sys.argv[2]
        # print node.name, node._count, \
        #     len(gt.get_leaves_by_name(node.name))

        # st.display_with_gene_tree(gt)

# st.display_with_gene_tree(gt)
gt._just_gains(st)
for node in st.traverse():
    print node._gains, node._justgains
    try:
        assert node._gains == node._justgains
    except:
        print 'AssertionError'
        break
