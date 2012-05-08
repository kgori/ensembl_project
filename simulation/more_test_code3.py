#!/usr/bin/env python

import ete2, re, sys
from ete2.parser.newick import read_newick, write_newick
from collections import defaultdict
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

    def _get_species_subtree(self,species_tree):
        """
        Recursive function for defining the subtree of the species tree associated
        with each node of the gene tree
        """
        # BASE CASE
        if not self.children:
            self._species_subtree = species_tree.get_leaves_by_name(self.name)[0]
            return
        # VARIOUS OPTIONS FOR SETTING SPECIES SUBTREE
        while True:
            if self._species_subtree: # subtree already set by being passed from parent speciation node
                break
            if self.type()=="S":
                self._species_subtree = species_tree.get_common_ancestor(self._taxon_set)
                for child in self.children:
                    if child.type()=="D":
                        for schild in self._species_subtree.children:
                            if child._taxon_set.issubset(schild._taxon_set):
                                child._species_subtree = schild
                break
            if not self.up: # root node
                self._species_subtree = species_tree.get_common_ancestor(self._taxon_set)
                break
            if self.type()=="D":
                self._species_subtree = self.up._species_subtree
                break 
        # RECURSIVE CASE
        for child in self.children:
            child._get_species_subtree(species_tree)

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
        self.add_features(name='NoName')
        self.features.update(['dist', 'support'])

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
            if n.type()=='D':
                nstyle['fgcolor']='red'
                nstyle['size']=8
            elif n.type()=='S':
                nstyle['fgcolor']='blue'
                nstyle['size']=8
            else: 
                nstyle['fgcolor']='lightblue'
            n.set_style(nstyle)
        self.show()


class SpeciesTree(ete2.Tree):

    def parse_species_name(self):
        if not self.children:
            try: return re.compile('(?<=_)[\w]+').search(self.name).group()
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
        self._already_deleted = []  # added

        self.features = set([])
        # Add basic features
        self.add_features(name='NoName')
        self.features.update(['dist', 'support'])

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
                if n._gains>0: n.add_face(ete2.TextFace(n._gains,fgcolor='green'),column=0,\
                    position='branch-top')
            if hasattr(n,'_losses'):
                if n._losses: n.add_face(ete2.TextFace(n._losses,fgcolor='red'),column=0,\
                    position='branch-bottom')
        self.show()

    def display_with_gene_tree(self,genetree,scores=True):
        for n in self.traverse():
            if hasattr(n,'_gains'):
                if n._gains>0: n.add_face(ete2.TextFace('{0:+}'.format(n._gains),fgcolor='green'),\
                    column=0,position='branch-top')
            if hasattr(n,'_losses'):
                if n._losses: n.add_face(ete2.TextFace('{0:+}'.format(-1*n._losses),fgcolor='red'),\
                    column=0,position='branch-bottom')
            if hasattr(n,'_count'):
                if n._count: n.add_face(ete2.TextFace('{0}'.format(n._count),fgcolor='blue'),\
                    column=0,position='branch-right')
            nstyle=ete2.NodeStyle()
            if n._events:
                if 'Earliest evidence for this gene from init' in n._events: nstyle['fgcolor']='cyan'
            n.set_style(nstyle)
        for n in genetree.traverse():
            nstyle = ete2.NodeStyle()
            if n.type()=='D':
                nstyle['fgcolor']='red'
                nstyle['size']=6
            elif n.type()=='S':
                nstyle['fgcolor']='blue'
                nstyle['size']=6
            else: 
                nstyle['fgcolor']='lightblue'
            n.set_style(nstyle)
       
        def layout(node):
            species_ts = ete2.TreeStyle()
            species_ts.scale=6
            species_ts.branch_vertical_margin=1.5
            if not node.up:
                ete2.faces.add_face_to_node(ete2.TreeFace(self,species_ts),node,0,position='branch-top')

        ts = ete2.TreeStyle()
        ts.layout_fn=layout
        ts.title.add_face(ete2.TextFace(sys.argv[2],fsize=20),column=0)
        genetree.show(tree_style=ts)

    def _find_losses(self, gene_tree_node, species_subtree=None, init=False,\
        already_deleted=[]):
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
            for node in deletions:
                node._events.append('Gene loss implied by init')
                node._losses += 1
        else:
            for node in deletions:
                node._events.append('Gene loss implied by duplication node')
                node._losses += 1

        return deletions + already_deleted

    def _add_duplication_nodes(self, gene_tree_node, already_deleted = [], verbose=False, **kwargs):
        # CASE 1 - root node
        if not gene_tree_node.up:
            if gene_tree_node.type()=='D':   
                if verbose: print 'DUPLICATION NODE FOUND'
                if verbose: print gene_tree_node, already_deleted
                mrca = self._find_mrca(gene_tree_node._taxon_set)
                mrca._gains += 1
                mrca._events.append('Gene duplication')
                for child in gene_tree_node.children:
                    self._add_duplication_nodes(child, already_deleted=already_deleted)
            else:
                if verbose: print 'ROOT NODE FOUND, RECURSING'
                for child in gene_tree_node.children:
                    self._add_duplication_nodes(child, already_deleted=already_deleted) 
        # CASE 2 - duplication following duplication                                    # This case checks for speciation followed by loss,
        elif gene_tree_node.type()=='D' and gene_tree_node.up.type()=='D':              # which would result in a missing speciation node
            if verbose: print 'DUPLICATION FOLLOWING DUPLICATION',already_deleted       # that would otherwise be observed between the
            mrca = self._find_mrca(gene_tree_node._taxon_set)                           # duplication nodes
            mrca._gains += 1
            mrca._events.append('Gene duplication') 
            # subtree = gene_tree_node._species_subtree 
            subtree = self._find_mrca(gene_tree_node.up._taxon_set)                  
            already_deleted = self._find_losses(gene_tree_node, subtree,\
                already_deleted=already_deleted)
            for d in already_deleted:
                if d not in self._already_deleted:
                    if verbose: print 'adding', d
                    self._already_deleted.append(d)           
            for child in gene_tree_node.children:
                self._add_duplication_nodes(child, already_deleted=already_deleted)        
        # CASE 3 - duplication                                                          # Record gain at mrca, then recurse on child nodes
        elif gene_tree_node.type()=='D':                                                
            if verbose: print 'DUPLICATION NODE FOUND'
            if verbose: print gene_tree_node, already_deleted
            mrca = self._find_mrca(gene_tree_node._taxon_set)
            mrca._gains += 1
            mrca._events.append('Gene duplication')
            for child in gene_tree_node.children:
                self._add_duplication_nodes(child, already_deleted=already_deleted)
        # CASE 4 - speciation following duplication                                     # Using mrca found at parent duplication node,
        elif gene_tree_node.type()!='D' and gene_tree_node.up.type()=='D':              # record deletions under this node in already_
            if verbose: print 'SPECIATION FOLLOWING DUPLICATION',already_deleted        #-deleted list and _already_deleted global,  
            subtree = self._find_mrca(gene_tree_node.up._taxon_set)                     # then recurse on child nodes
            already_deleted = self._find_losses(gene_tree_node, subtree,\
                already_deleted=already_deleted)
            for d in already_deleted:
                if d not in self._already_deleted:
                    if verbose: print 'adding', d
                    self._already_deleted.append(d)
            if gene_tree_node.children:
                for child in gene_tree_node.children:
                    self._add_duplication_nodes(child, already_deleted=already_deleted)
            else: return
        # CASE 5 - leaf node                                                            # base case, end recursion
        elif not gene_tree_node.children:                                               
            if verbose: print 'LEAF NODE FOUND, RETURNING'
            if verbose: print gene_tree_node, already_deleted
            return   
        # CASE 6 - anything else
        else:                                                                           # Just recurse on child nodes
            if verbose: print 'SPECIATION NODE FOUND', already_deleted
            if verbose: print gene_tree_node
            for child in gene_tree_node.children:
                self._add_duplication_nodes(child, already_deleted=already_deleted)
       
    def _add_root_node(self, gene_tree_node):
        # if not gene_tree_node.type() == 'D':
        # print 'ADDING ROOT NODE'
        mrca = self._find_mrca(gene_tree_node._taxon_set)
        mrca._count += 1
        mrca._events.append('Earliest evidence for this gene from init')
        self._find_losses(gene_tree_node, init=True, already_deleted=self._already_deleted)

    def _finalise_counts(self):
        self._count += (self._gains - self._losses)
        for node in self.iter_descendants('preorder'):
            node._count += (node.up._count + node._gains - node._losses)

    def add_gene_tree(self, gene_tree):
        # gene_tree._get_species_subtree(self)
        self._add_duplication_nodes(gene_tree)
        self._add_root_node(gene_tree)
        self._finalise_counts()
        self._already_deleted = []


species_tree_newick=open(sys.argv[1]).read()
gene_tree_newick = open(sys.argv[2]).read()
sp_tr = SpeciesTree(species_tree_newick,format=1)
gn_tr = GeneTree(gene_tree_newick,format=1)
gn_tr.get_descendant_evol_events()
sp_tr.add_gene_tree(gn_tr)
for taxon in sp_tr._taxon_set:
    print taxon, (sp_tr&taxon)._count, len(gn_tr.get_leaves_by_name(taxon))#len([n.name.split('_')[-1] for n in gn_tr.get_leaves() if n.name.split('_')[-1]==taxon])
    # if (sp_tr&taxon)._count != len([n.name.split('_')[-1] for n in gn_tr.get_leaves() if n.name.split('_')[-1]==taxon]): print taxon
    assert (sp_tr&taxon)._count == len(gn_tr.get_leaves_by_name(taxon))#len([n.name.split('_')[-1] for n in gn_tr.get_leaves() if n.name.split('_')[-1]==taxon])
# sp_tr.display_with_gene_tree(gn_tr)
# print sp_tr._already_deleted
# st = SpeciesTree('(SE001:60.971176,(SE002:54.097175,((SE003:8.863838,SE013:8.7271048):33.859374,(((SE004:40.303116,(SE007:35.957807,(SE008:27.55188,(SE010:23.561527,SE011:23.148886):2.2511886):0.10991615):3.753085):7.938699,SE006:35.196556):7.992514,((SE005:16.104527,(SE014:1.9441679,SE015:3.2796744):19.411262):15.863363,(SE009:19.907848,SE012:39.749202):5.4071946):19.687332):3.7807204):17.137093):7.2302568):0;')
# gt=GeneTree('(SE001,(SE002,((((SE004,(SE007,(SE008,SE010))),SE006),((SE004,(SE007,(SE008,(SE008,SE010))[&&NHX:D=Y])),SE006))[&&NHX:D=Y],(((SE005,(SE014,SE015)),(SE009,SE012)),(SE005,(SE009,SE012)))[&&NHX:D=Y])));',format=1)
# st2 =SpeciesTree('(((1,2),(3,(4,5))),(6,7));')
# gt2 =GeneTree('((((1,2),(4,5)),((1,2),(4,5)))[&&NHX:D=Y](6,7));')
# st2.add_gene_tree(gt2)
# print st2._already_deleted
# st._initialise(gt)
# st.add_gene_tree(gt)
# st2.display_with_gene_tree(gt2)
# st.display_with_gene_tree(gt)
# big_st = "(((Hsap,Ptro),Nleu),(Mmus,Rnor));"
# big_gt = "(((((gene1_Hsap,gene1_Ptro)[&&NHX:D=N],gene1_Nleu)[&&NHX:D=N],((gene1a_Mmus,gene1a_Rnor)[&&NHX:D=N],(gene1b_Mmus,gene1b_Rnor)[&&NHX:D=N])[&&NHX:D=Y])[&&NHX:D=N],((((gene2a_Hsap,gene2a_Ptro)[&&NHX:D=N],gene2a_Nleu)[&&NHX:D=N],((gene2b_Hsap,gene2c_Hsap)[&&NHX:D=Y],gene2b_Ptro)[&&NHX:D=N])[&&NHX:D=Y],((gene2a_Mmus,gene2a_Rnor)[&&NHX:D=N],(gene2b_Mmus,gene2b_Rnor)[&&NHX:D=N])[&&NHX:D=Y])[&&NHX:D=N])[&&NHX:D=Y],(((gene3_Hsap,gene3_Ptro)[&&NHX:D=N],gene3_Nleu)[&&NHX:D=N],((gene3a_Mmus,gene3_Rnor)[&&NHX:D=N],gene3b_Mmus)[&&NHX:D=Y])[&&NHX:D=N])[&&NHX:D=Y];"
# big_st = SpeciesTree(big_st)
# big_gt = GeneTree(big_gt)
# big_st.add_gene_tree(big_gt)
# big_st.display_with_gene_tree(big_gt)
# small_st = "((Hsap,Ptro),Mmus);"
# small_gt = "(((gene1_Hsap,gene1_Mmus)[&&NHX:D=N],(((gene2a_Hsap,gene2a_Ptro)[&&NHX:D=N],(gene2b_Hsap,gene2b_Ptro)[&&NHX:D=N])[&&NHX:D=Y],(gene2a_Mmus,gene2b_Mmus)[&&NHX:D=Y])[&&NHX:D=N])[&&NHX:D=Y],((gene3_Hsap,gene3_Ptro)[&&NHX:D=N],(gene3_Mmus,gene3'_Mmus)[&&NHX:D=Y])[&&NHX:D=N])[&&NHX:D=Y];"
# compara_st = "(Hsap,Mmus);"
# compara_gt = "(((gene1_Hsap,gene1_Mmus)[&&NHX:D=N],((gene2a_Hsap,gene2b_Hsap)[&&NHX:D=Y],(gene2a_Mmus,gene2b_Mmus)[&&NHX:D=Y])[&&NHX:D=N])[&&NHX:D=Y],(gene3_Hsap,(gene3a_Mmus,gene3b_Mmus)[&&NHX:D=Y])[&&NHX:D=N])[&&NHX:D=Y];"
# small_st = SpeciesTree(small_st)
# small_gt = GeneTree(small_gt)
# compara_gt = GeneTree(compara_gt)
# compara_st = SpeciesTree(compara_st)
# compara_st.add_gene_tree(compara_gt)
# small_st.add_gene_tree(small_gt)
# small_st.display_with_gene_tree(small_gt)
# compara_st.display_with_gene_tree(compara_gt)
