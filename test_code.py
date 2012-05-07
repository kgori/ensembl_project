#!/usr/bin/env python

import ete2, re
from ete2.parser.newick import read_newick, write_newick
try:
    from ete2.treeview.main import NodeStyle, _FaceAreas, FaceContainer, FACE_POSITIONS
    from ete2.treeview.faces import Face
except ImportError:
    TREEVIEW = False
else:
    TREEVIEW = True


class Gene_Tree(ete2.Tree):

    def __init__(self, newick=None, format=0):
        ### From ete2.Tree class definition
        self._children = []
        self._up = None
        self._dist = 1.0
        self._support = 1.0
        self._rootdist = 0 #added

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
            n=node
            while node.up:
                if node.up.type()=='N':
                    n._rootdist +=1
                node=node.up

    def parse_species_name(self):
        if not self.children:
            return re.compile("(?<=_)[\w]+").search(self.name).group()

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
        if hasattr(self,'D'): return self.D
        else: return 'L'

    def get_species(self):
        for x in self.iter_leaves():
            yield x.parse_species_name()

    def display(self):
        for n in self.traverse():
            print n.type(),
            nstyle = ete2.NodeStyle()
            if n.type()=="Y":
                nstyle["fgcolor"]="red"
                nstyle["size"]=8
                print "red"
            elif n.type()=="N":
                nstyle["fgcolor"]="blue"
                nstyle["size"]=8
                print "blue"
            else: 
                nstyle["fgcolor"]="lightblue"
                print "lightblue"
            n.set_style(nstyle)
        self.show()


class Species_Tree(ete2.Tree):

    def __init__(self, newick=None, format=0):
        ### From ete2.Tree class definition
        self._children = []
        self._up = None
        self._dist = 1.0
        self._support = 1.0
        self._count = 0             # added
        self._change = 0            # added
        self._rootdist = 0          # added
        self._fitchset = set([])    # added

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
            n=node
            while node.up:
                n._rootdist +=1
                node=node.up

    def display(self):
        for n in self.traverse():
            if hasattr(n,'count'):
                n.add_face(ete2.TextFace(n.count),column=0,position="branch-top")
        self.show()
    
    def reset(self):
        for n in self.traverse():
            n._count=0
            n._change=0
            n._fitchset=set([])

    def _render(self,filename,**kwargs):
        self.render(filename,**kwargs)

    def display_with_gene_tree(self,genetree,scores=True):
        if scores:
            for n in self.traverse():
                n.add_face(ete2.TextFace(n._count),column=0,position="branch-top")
                n.add_face(ete2.TextFace("{0:+}".format(n._change)),column=1,position="branch-bottom")
        for n in genetree.traverse():
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
       
        def layout(node):
            species_ts = ete2.TreeStyle()
            species_ts.scale=60
            species_ts.branch_vertical_margin=15
            if not node.up:
                ete2.faces.add_face_to_node(ete2.TreeFace(self,species_ts),node,0,position="branch-top")

        ts = ete2.TreeStyle()
        ts.layout_fn=layout
        genetree.show(tree_style=ts)

    def render_with_gene_tree(self,genetree,filename,scores=True, **kwargs):
        if scores:
            for n in self.traverse():
                n.add_face(ete2.TextFace(n._count),column=0,position="branch-top")
                n.add_face(ete2.TextFace("{0:+}".format(n._change)),column=1,position="branch-bottom")
        for n in genetree.traverse():
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

        
        def layout(node):
            species_ts = ete2.TreeStyle()
            species_ts.scale=60
            if not node.up:
                ete2.faces.add_face_to_node(ete2.TreeFace(self,species_ts),node,0,position="branch-top")

        ts = ete2.TreeStyle()
        ts.layout_fn=layout
        genetree.render(filename,tree_style=ts, **kwargs)

    def _add(self,amount):
        self._count+=amount

    def _add_counts(self, species_set, amount=1):
        self._run_fitch_algorithm(species_set)
        self._add_fitch_scores(amount)

    def _run_fitch_algorithm(self,species_set):
        """
        Given an input node, determine the set of states available at the input node, which is the intersection
        of the states of its child nodes, or the union if the intersection is empty.
        Start at the tips and work back.
        """
        for node in self.traverse('postorder'):
            if node.children:
                sets = [child._fitchset for child in node.children]
                i = self._multiintersection(sets)
                if i:
                    node._fitchset = i
                else: node._fitchset = self._multiunion(sets)
            else:
                if (set([node.name]) & species_set):
                    node._fitchset=set([1])
                else: node._fitchset=set([0])

    def _add_fitch_scores(self,amount=1):
        for node in self.traverse('postorder'):
            try:
                if node._fitchset == set([1]):
                    node._add(amount)
                elif node._fitchset == set([0,1]):
                    if node.up:
                        if node.up._fitchset == set([1]):
                            node._add(amount)
                node._fitchset = set([])
            except AttributeError: # Trying to access a fitchset too deep in the tree, where it is None
                continue

    def _initialise(self,gene_tree):
        """
        Adds count of 1 to all nodes in species tree which lead to a present species in the
        gene tree
        """
        gene_tree_species_set = set(gene_tree.get_species())
        self._run_fitch_algorithm(gene_tree_species_set)
        self._add_fitch_scores(1)

    def _multiunion(self,l):
        """ 
        Returns the union of all sets in a given list of sets
        Adapted from Learning Python, Mark Lutz, O'Reilly Media
        """
        res = []
        for s in l:
            for x in s:
                if x not in res:
                    res.append(x)
        return set(res)

    def _multiintersection(self,l):
        """ 
        Returns the intersection of all sets in a list
        """
        i = l[0]
        for s in l[1:]:
            i = i & s
        return i

    def add_gene_tree(self,gene_tree):
        self._initialise(gene_tree)
        for node in gene_tree.traverse():
            if node.type() == 'Y':
                sets = [set(child.get_species()) for child in node.children]
                union = self._multiunion(sets)     

                for s,species_subtree in zip(sets,[self.get_common_ancestor(s) if len(s)>1 else self.get_leaves_by_name(list(s)[0])[0] for s in sets ]):
                    number_of_speciations = species_subtree._rootdist - node._rootdist
                    for i in range(number_of_speciations):
                        species_subtree=species_subtree.up
                    species_subtree._add_counts(s,1)
                               
                species_subtree._add_counts(union,-1)
        for node in self.traverse():
            if node.up: node._change=node._count-node.up._count
            else: node._change=node._count
        
    def print_count(self):
        for node in self.traverse('postorder'):
            print node,node._count



small_st = "((Hsap,Ptro),Mmus);"
small_gt = "(((gene1_Hsap,gene1_Mmus)[&&NHX:D=N],(((gene2a_Hsap,gene2a_Ptro)[&&NHX:D=N],(gene2b_Hsap,gene2b_Ptro)[&&NHX:D=N])[&&NHX:D=Y],(gene2a_Mmus,gene2b_Mmus)[&&NHX:D=Y])[&&NHX:D=N])[&&NHX:D=Y],((gene3_Hsap,gene3_Ptro)[&&NHX:D=N],(gene3_Mmus,gene3'_Mmus)[&&NHX:D=Y])[&&NHX:D=N])[&&NHX:D=Y];"
compara_st = "(Hsap,Mmus);"
compara_gt = "(((gene1_Hsap,gene1_Mmus)[&&NHX:D=N],((gene2a_Hsap,gene2b_Hsap)[&&NHX:D=Y],(gene2a_Mmus,gene2b_Mmus)[&&NHX:D=Y])[&&NHX:D=N])[&&NHX:D=Y],(gene3_Hsap,(gene3a_Mmus,gene3b_Mmus)[&&NHX:D=Y])[&&NHX:D=N])[&&NHX:D=Y];"
big_st = "(((Hsap,Ptro),Nleu),(Mmus,Rnor));"
big_gt = "(((((gene1_Hsap,gene1_Ptro)[&&NHX:D=N],gene1_Nleu)[&&NHX:D=N],((gene1a_Mmus,gene1a_Rnor)[&&NHX:D=N],(gene1b_Mmus,gene1b_Rnor)[&&NHX:D=N])[&&NHX:D=Y])[&&NHX:D=N],((((gene2a_Hsap,gene2a_Ptro)[&&NHX:D=N],gene2a_Nleu)[&&NHX:D=N],((gene2b_Hsap,gene2c_Hsap)[&&NHX:D=Y],gene2b_Ptro)[&&NHX:D=N])[&&NHX:D=Y],((gene2a_Mmus,gene2a_Rnor)[&&NHX:D=N],(gene2b_Mmus,gene2b_Rnor)[&&NHX:D=N])[&&NHX:D=Y])[&&NHX:D=N])[&&NHX:D=Y],(((gene3_Hsap,gene3_Ptro)[&&NHX:D=N],gene3_Nleu)[&&NHX:D=N],((gene3a_Mmus,gene3_Rnor)[&&NHX:D=N],gene3b_Mmus)[&&NHX:D=Y])[&&NHX:D=N])[&&NHX:D=Y];"
big_st = Species_Tree(big_st)
big_gt = Gene_Tree(big_gt)
small_st = Species_Tree(small_st)
small_gt = Gene_Tree(small_gt)
compara_gt = Gene_Tree(compara_gt)
compara_st = Species_Tree(compara_st)
compara_st.add_gene_tree(compara_gt)
small_st.add_gene_tree(small_gt)
compara_st.display_with_gene_tree(compara_gt,False)
compara_st.display_with_gene_tree(compara_gt,True)
big_st.add_gene_tree(big_gt)
# big_st.print_count()
# big_st.display()
# small_st.display()
# compara_st.display()
# big_gt.display()
# small_gt.display()
# compara_gt.display()
# big_st.render_with_gene_tree(big_gt,"big_both.pdf")
# small_st.render_with_gene_tree(small_gt,"small_both.pdf")
# compara_st.render_with_gene_tree(compara_gt,"compara_both.pdf")
