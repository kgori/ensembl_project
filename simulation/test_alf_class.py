#!/usr/bin/env python

import ete2, re, sys
from ete2.parser.newick import read_newick, write_newick
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
        self._taxon_set = None          # added
        self._species_subtree = None    # added
        self._index = None              # added

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
        
    def _render(self,filename,**kwargs):
        for n in self.traverse():
            nstyle = ete2.NodeStyle()
            if n.type()=="Y":
                nstyle["fgcolor"]="red"
                nstyle["size"]=8
            elif n.type()=="N":
                nstyle["fgcolor"]="blue"
                nstyle["size"]=8
            else: 
                nstyle["fgcolor"]="lightblue"
            n.set_style(nstyle)

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
                   
    def iter_duplication_nodes(self):
        for node in self.traverse("levelorder"):
            if node.type()=="D":
                yield node 

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
        self._change = 0            # added
        self._fitchset = set([])    # added
        self._sankoff = [None,None] # added
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

    def display(self):
        for n in self.traverse():
            if hasattr(n,'count'):
                n.add_face(ete2.TextFace(n.count),column=0,position="branch-top")
        self.show()

    def _render(self,filename,**kwargs):
        """ """
        self.render(filename,**kwargs)

    def display_with_gene_tree(self,genetree,scores=True):
        if scores:
            for n in self.traverse():
                n.add_face(ete2.TextFace(n._count),column=0,position="branch-top")
                n.add_face(ete2.TextFace("{0:+}".format(n._change)),column=1,position="branch-bottom")
                nstyle=ete2.NodeStyle()
                if n._events:
                    if "Initialised" in n._events: nstyle["fgcolor"]="cyan"
                    else: nstyle["fgcolor"]="green"
                    nstyle["size"]=6*len(n._events)
                elif n._change < 0:
                    nstyle["fgcolor"]="yellow"
                    nstyle["size"]=-6*n._change
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

    def render_with_gene_tree(self,genetree,filename,scores=True, **kwargs):
        if scores:
            for n in self.traverse():
                n.add_face(ete2.TextFace(n._count),column=0,position="branch-top")
                n.add_face(ete2.TextFace("{0:+}".format(n._change)),column=1,position="branch-bottom")
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
            species_ts.scale=60
            if not node.up:
                ete2.faces.add_face_to_node(ete2.TreeFace(self,species_ts),node,0,position="branch-top")

        ts = ete2.TreeStyle()
        ts.layout_fn=layout
        genetree.render(filename,tree_style=ts, **kwargs)

    def _add(self,amount):
        """ """
        self._count+=amount

    def _run_sankoff_algorithm(self,amount,taxon_set,c=[[0,2],[1,0]]):
        # BACKWARDS PASS
        for node in self.traverse('postorder'):
            if not node.children:
                if (set([node.name]) & taxon_set):
                    node._sankoff=[sys.maxint,0]
                else: 
                    node._sankoff = [0,sys.maxint]
            else:
                node._sankoff = [min([c[i][j] + node.children[0]._sankoff[j] for j in range(2)]) +\
                                    min([c[i][k] + node.children[1]._sankoff[k] for k in range(2)])\
                                    for i in range(2)]

        # FORWARDS PASS
        for node in self.traverse('preorder'):
            if node._sankoff[1] < node._sankoff[0]:
                node._add(amount)
                
        return min(self._sankoff)

    def _run_fitch_algorithm(self,amount,taxon_set):
        """
        Given an input node, determine the set of states available at the input node, which is the intersection
        of the states of its child nodes, or the union if the intersection is empty.
        Start at the tips and work back.
        """
        score = 0
        # BACKWARDS PASS
        for node in self.traverse('postorder'):
            if node.children:
                sets = [child._fitchset for child in node.children]
                i = self._multiintersection(sets)
                if i:
                    node._fitchset = i
                else: 
                    node._fitchset = self._multiunion(sets)
                    score += 1
            else:
                if (set([node.name]) & taxon_set):
                    node._fitchset=set([1])
                else: node._fitchset=set([0])

        # FORWARDS PASS
        for node in self.traverse('preorder'):
            if node._fitchset == set([0,1]):
                #root case
                if node == self:
                    node._fitchset = set([1]) # ACCTRAN kind of behaviour
                    node._add(amount)
                else:
                    if node.up._fitchset == set([1]):
                        node._fitchset = set([1])
                        node._add(amount)
                    elif node.up._fitchset == set([0]):
                        node._fitchset = set([0])
                    
            elif node._fitchset == set([1]):
                node._add(amount)
        
        return score
            
    def _initialise(self,gene_tree):
        """
        Adds count of 1 to all nodes in species tree which lead to a present species in the
        gene tree
        """
        mrca = self._find_mrca(gene_tree._taxon_set)
        mrca._events.append("Initialised")
        # mrca._run_fitch_algorithm(1,gene_tree._taxon_set)
        self._run_sankoff_algorithm(1,gene_tree._taxon_set)

    def _multiunion(self,l):
        """ 
        Returns the union of all sets in a given list of sets
        """
        u = l[0]
        for s in l[1:]:
            u = u | s
        return u

    def _multiintersection(self,l):
        """ 
        Returns the intersection of all sets in a given list of sets
        """
        i = l[0]
        for s in l[1:]:
            i = i & s
        return i                 

    def _find_mrca(self,taxon_set):
        """
        Returns the most recent common ancestor of the taxon set of the node
        """
        if len(taxon_set)>1:
            return self.get_common_ancestor(taxon_set)
        else:
            return self.get_leaves_by_name(list(taxon_set)[0])[0]

    def add_gene_tree(self,gene_tree):
        # gene_tree._get_species_subtree(self)
        self._initialise(gene_tree)
        index=0
        for node in gene_tree.iter_duplication_nodes():
            index+=1
            mrca = self._find_mrca(node._taxon_set)
            mrca._events.append("Gain here as mrca of {0}".format(index))
            mrca._run_fitch_algorithm(-1,node._taxon_set)
            for child in node.children:
                mrca._run_fitch_algorithm(1,child._taxon_set)
                
        for node in self.traverse():
            if node.up: node._change=node._count-node.up._count
            else: node._change=node._count



species_tree_newick=open(sys.argv[1]).read()
sp_tr = SpeciesTree(species_tree_newick,format=1)
gene_tree_newick = open(sys.argv[2]).read()
gn_tr = GeneTree(gene_tree_newick,format=1)
gn_tr.get_descendant_evol_events()
sp_tr.add_gene_tree(gn_tr)

for node in sp_tr.traverse():
    print node,node._sankoff,node._count,node._events

sp_tr.display_with_gene_tree(gn_tr)



# OBSOLESCENCE

# def _get_species_subtree(self,species_tree):
#         """
#         Recursive function for defining the subtree of the species tree associated
#         with each node of the gene tree
#         """
#         # BASE CASE
#         if not self.children:
#             self._species_subtree = species_tree.get_leaves_by_name(self.name)[0]
#             return
#         # VARIOUS OPTIONS FOR SETTING SPECIES SUBTREE
#         while True:
#             if self._species_subtree: # subtree already set by being passed from parent speciation node
#                 break
#             if self.type()=="S":
#                 self._species_subtree = species_tree.get_common_ancestor(self._taxon_set)
#                 for child in self.children:
#                     if child.type()=="D":
#                         for schild in self._species_subtree.children:
#                             if child._taxon_set.issubset(schild._taxon_set):
#                                 child._species_subtree = schild
#                 break
#             if not self.up: # root node
#                 self._species_subtree = species_tree.get_common_ancestor(self._taxon_set)
#                 break
#             if self.type()=="D":
#                 self._species_subtree = self.up._species_subtree
#                 break 
#         # RECURSIVE CASE
#         for child in self.children:
#             child._get_species_subtree(species_tree)

# def _apply_counts_backwards(self,species_set,amount=1):
#     """
#     Starting at the leaves, add 'amount' to the leaf node.
#     Then add 'amount' to parent nodes if at least one
#     child node has a count > 0.
#     If all children of a node carry 0 members of a gene family,
#     that node carries 0 members.
#     """
#     for node in self.traverse('postorder'):
#         if not node.children:
#             if set([node.name]) & species_set:
#                 node._add(amount)
#                 node._tracker=1

#         else:        
#             if sum([1 if child._tracker!=0 else 0 for child in node.children]): # evaluates to true iff at least one child has gained score > 0
#                 node._add(amount)
#                 node._tracker=1
#                 for child in node.children: 
#                     child._tracker=0

# def _add_fitch_scores(self,amount=1):
#     for node in self.traverse('postorder'):
#         try:
#             if node._fitchset == set([1]):
#                 node._add(amount)
#             elif node._fitchset == set([0,1]):
#                 if node.up:
#                     if node.up._fitchset == set([1]):
#                         node._add(amount)
#             node._fitchset = set([])
#         except AttributeError: # Trying to access a fitchset too deep in the tree, where it is None
#             continue

# def _add_counts(self, species_set, amount=1):
#     self._run_fitch_algorithm(species_set)
#     self._add_fitch_scores(amount)
