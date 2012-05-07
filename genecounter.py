#!/usr/bin/env python

import ete2, re, sys
from ete2.parser.newick import read_newick, write_newick
from collections import defaultdict, deque
try:
    from ete2.treeview.main import NodeStyle, _FaceAreas, FaceContainer, FACE_POSITIONS
    from ete2.treeview.faces import Face
except ImportError:
    TREEVIEW = False
else:
    TREEVIEW = True

class GeneTree(ete2.PhyloTree):

    def parse_species_name(self):
        if not self.children:
            try: return self.name.split('_')[-1]
            except AttributeError: return self.name

    def get_species(self):
        for x in self.iter_leaves():
            yield x.parse_species_name()

    def __init__(self, newick=None, format=0):
        ### From ete2.Tree class definition
        self._children = []
        self._up = None
        self._dist = 1.0
        self._support = 1.0
        self._taxon_set = None              # added
        self._species_subtree = None        # added
        self._index = None                  # added

        self.features = set([])
        # Add basic features
        self.add_features(name="NoName")
        self.features.update(["dist", "support"])

        # Initialize tree
        if newick is not None:
            read_newick(newick, root_node = self, format=format)
            
        if TREEVIEW:
            self._faces = _FaceAreas()
        for node in self.traverse():
            node._taxon_set = set(node.get_species())

        def nf(name):
            return name
        self.set_species_naming_function(nf)

    def type(self):
        try:
            if hasattr(self,'evoltype'): return self.evoltype
            elif hasattr(self,'D'):
                if self.children:
                    if self.D=='Y': return 'D'
                    elif self.D=='N': return 'S'
                else: return 'L'
            else: return 'L'
        except AttributeError:
            return 'L'

    def display(self):
        for n in self.traverse():
            print n.type(),
            nstyle = ete2.NodeStyle()
            if n.type()=="D":
                nstyle["fgcolor"]="red"
                nstyle["size"]=8
            elif n.type()=="S":
                nstyle["fgcolor"]="blue"
                nstyle["size"]=8
            else: 
                nstyle["fgcolor"]="lightblue"
            n.set_style(nstyle)
        self.show()


class SpeciesTree(ete2.Tree):

    def parse_species_name(self):
        if not self.children:
            try: return re.compile("(?<=_)[\w]+").search(self.name).group()
            except AttributeError: return self.name
    
    def get_species(self):
        for x in self.iter_leaves():
            yield x.parse_species_name()

    def __init__(self, newick=None, format=0):
        ### From ete2.Tree class definition
        self._children = []
        self._up = None
        self._dist = 1.0
        self._support = 1.0
        self._count = 0             # added
        self._gains = 0             # added
        self._losses = 0            # added
        self._change = 0            # added
        self._events = []           # added
        self._taxon_set = None      # added

        self.features = set([])
        # Add basic features
        self.add_features(name="NoName")
        self.features.update(["dist", "support"])

        # Initialize tree
        if newick is not None:
            read_newick(newick, root_node = self, format=format)
            
        if TREEVIEW:
            self._faces = _FaceAreas()
        for node in self.traverse():
            node._taxon_set = set(node.get_species())

    def traceback(self, node):
        while node.up and self != node:
            yield node
            node = node.up
        yield node  

    def _find_mrca(self,taxon_set):
        """
        Returns the most recent common ancestor of the taxon set of the node
        """
        if len(taxon_set)>1:
            return self.get_common_ancestor(taxon_set)
        else:
            return self.get_leaves_by_name(list(taxon_set)[0])[0]  

    def display(self):
        for n in self.traverse():
            if hasattr(n,'_gains'):
                if n._gains>0: n.add_face(ete2.TextFace(n._gains,fgcolor="green"),column=0,position="branch-top")
            if hasattr(n,'_losses'):
                if n._losses: n.add_face(ete2.TextFace(n._losses,fgcolor="red"),column=0,position="branch-bottom")
        self.show()


    def display_with_gene_tree(self,genetree,scores=True):
        for n in self.traverse():
            if hasattr(n,'_gains'):
                if n._gains>0: n.add_face(ete2.TextFace('{0:+}'.format(n._gains),fgcolor="green"),column=0,position="branch-top")
            if hasattr(n,'_losses'):
                if n._losses: n.add_face(ete2.TextFace('{0:+}'.format(-1*n._losses),fgcolor="red"),column=0,position="branch-bottom")
            if hasattr(n,'_count'):
                if n._count: n.add_face(ete2.TextFace('{0}'.format(n._count),fgcolor="blue"),column=0,position="branch-right")
            nstyle=ete2.NodeStyle()
            if n._events:
                if "Earliest evidence for this gene from init" in n._events: nstyle["fgcolor"]="cyan"
                else: nstyle["fgcolor"]="green"
                nstyle["size"]=6
            elif n._change < 0:
                nstyle["fgcolor"]="yellow"
                nstyle["size"]=6
            n.set_style(nstyle)
        for n in genetree.traverse():
            nstyle = ete2.NodeStyle()
            if n.type()=="D":
                nstyle["fgcolor"]="red"
                nstyle["size"]=8
            elif n.type()=="S":
                nstyle["fgcolor"]="blue"
                nstyle["size"]=8
            else: 
                nstyle["fgcolor"]="lightblue"
            n.set_style(nstyle)
       
        def layout(node):
            species_ts = ete2.TreeStyle()
            species_ts.scale=6
            species_ts.branch_vertical_margin=1.5
            if not node.up:
                ete2.faces.add_face_to_node(ete2.TreeFace(self,species_ts),node,0,position="branch-top")

        ts = ete2.TreeStyle()
        ts.layout_fn=layout
        ts.title.add_face(ete2.TextFace(sys.argv[2],fsize=20),column=0)
        genetree.show(tree_style=ts)

    def _initialise(self, gene_tree_node, **kwargs):
        if not gene_tree_node.type() == "D":
            return self.find_events(gene_tree_node, **kwargs)
        else:
            return set([])

    def find_events(self, gene_tree_node, species_subtree=None, init=False, gains=True, losses=True, already_deleted=[]):
        """
        inset = set of leaf nodes present in the genetree
        """
        inset = gene_tree_node._taxon_set
        if not species_subtree: species_subtree = self._find_mrca(inset)
        totalset = set(species_subtree.get_leaf_names())
        outset = totalset-inset
        traced = []
        untraced = defaultdict(list)
        deletions = []
        
        for inleaf in inset:
            inNode = self&inleaf
            for n in species_subtree.traceback(inNode):
                if n not in traced:
                    traced.append(n)
                else: break

        for outleaf in outset:
            outNode = self&outleaf
            for n in species_subtree.traceback(outNode):
                if n not in traced:
                    untraced[outleaf].append(n)
                else: break        

        # each value in untraced is a list of all the untraced nodes
        # under each leaf in the outset. We only want to count the
        # earliest ones.
        for val in untraced.values(): 
            earliest = sorted(val, key =\
                lambda x: self.get_distance(x, topology_only = True))[0]
            if earliest not in (deletions + already_deleted):
                deletions.append(earliest)
        
        if init:
            if gains:
                species_subtree._events.append("Earliest evidence for this gene from init")
                species_subtree._count += 1
            if losses:
                for node in deletions:
                    node._events.append("Gene loss implied by init")
                    node._losses += 1
        else:
            if gains:
                species_subtree._events.append("Gene gain implied by duplication node")
                species_subtree._gains += 1
            if losses:
                for node in deletions:
                    node._events.append("Gene loss implied by duplication node")
                    node._losses += 1

        return deletions + already_deleted

    def add_gene_tree(self, gene_tree_node, already_deleted = [], **kwargs):
        # CASE 1 - leaf node
        if not gene_tree_node.children:  # base case,
            print "LEAF NODE FOUND, RETURNING"
            return                  # end recursion
        # CASE 2 - root node, either speciation or duplication
        elif not gene_tree_node.up:
            if not gene_tree_node.type()=="D":          
                print "INITIALISING ON ROOT NODE, RECURSING"                                       # Do initialisation stuff,
                already_deleted = self._initialise(gene_tree_node, init=True)   # then recurse on child nodes
                for child in gene_tree_node.children:
                    self.add_gene_tree(child, already_deleted=already_deleted)     
            else:
                print "DUPLICATION NODE FOUND AT ROOT"
                subtree = self._find_mrca(gene_tree_node._taxon_set)
                print "INITIALISING ON ROOT NODE"
                subtree._count += 1
                print "IMPUTING GAINS"
                already_deleted = self.find_events(gene_tree_node, losses=False, already_deleted=already_deleted)
                for child in gene_tree_node.children:
                    print "IMPUTING LOSSES, RECURSING"
                    already_deleted = self.find_events(child, species_subtree=subtree, gains=False, losses=True, already_deleted=already_deleted )
                    self.add_gene_tree(child, already_deleted=already_deleted)  
        # CASE 3 - duplication
        elif gene_tree_node.type()=="D":  # Do duplication node stuff, then recurse on child nodes
            print "DUPLICATION NODE FOUND"
            print gene_tree_node
            subtree = self._find_mrca(gene_tree_node._taxon_set)
            print "IMPUTING GAINS"
            already_deleted = self.find_events(gene_tree_node, losses=False, already_deleted=already_deleted)
            for child in gene_tree_node.children:
                print "IMPUTING LOSSES, RECURSING"
                already_deleted = self.find_events(child, species_subtree=subtree, gains=False, losses=True, already_deleted=already_deleted )
                self.add_gene_tree(child, already_deleted=already_deleted)
        # CASE 4 - speciation (or other)
        else:           # Just recurse on child nodes
            print "SPECIATION NODE FOUND, RECURSING"
            for child in gene_tree_node.children:
                self.add_gene_tree(child, already_deleted=already_deleted)
        
    def finalise_counts(self):
        for node in self.iter_descendants('preorder'):
            node._count += (node.up._count + node._gains - node._losses)