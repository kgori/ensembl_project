#!/usr/bin/env python

import ete2, re

class Gene_Tree(ete2.PhyloTree):
    def parse_species_name(self):
        if not self.children:
            return re.compile("(?<=_)[\w]+(?=_)").search(self.name).group()
        
    def type(self):
        if hasattr(self,'t'): return self.t
        else: return 'L'

    def get_species(self):
        for x in self.iter_leaves():
            yield x.parse_species_name()


class Species_Tree(ete2.PhyloTree):
    count=0 # Initial state is for each node in the species tree to have 0 genes

    def _add(self,amount):
        self.count+=amount

    def _apply_counts_backwards(self,species_set,amount=1):
        """
        Starting at the leaves, add 'amount' to the leaf node.
        Then add 'amount' to parent nodes if at least one
        child node has a count > 0.
        If all children of a node carry 0 members of a gene family,
        that node carries 0 members.
        """
        for node in self.traverse('postorder'):
            if node.children:
                if sum([1 if child.count>0 else 0 for child in node.children]): # evaluates to true iff at least one child has score > 0
                    node._add(amount)
            else:
                if set([node.name]) & species_set:
                    node._add(amount)

    def _initialise(self,gene_tree):
        """
        Adds count of 1 to all nodes in species tree which lead to a present species in the
        gene tree
        """
        gene_tree_species_set = set(gene_tree.get_species())
        self._apply_counts_backwards(gene_tree_species_set,1)

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

    def add_gene_tree(self,gene_tree):
        self._initialise(gene_tree)
        for node in gene_tree.traverse():
            if node.type() == 'D':
                sets = [set(child.get_species()) for child in node.children]
                union = self._multiunion(sets)     
                for s,t in zip(sets,[self.get_common_ancestor(s) if len(s)>1 else self.get_leaves_by_name(list(s)[0])[0] for s in sets ]):
                    t._apply_counts_backwards(s,1)
                if len(union)>1:
                    union_tree = self.get_common_ancestor(union)
                else: 
                    union_tree = self.get_leaves_by_name(list(union)[0])[0]
                union_tree._apply_counts_backwards(union,-1)

    def print_count(self):
        for node in self.traverse('postorder'):
            print node,node.count


species_tree = "(Hsap,Mmus);"
duplicates_tree = "(((1_Hsap_,1_Mmus_)[&&NHX:t=S],((2_Hsap_,2'_Hsap_)[&&NHX:t=D],(2_Mmus_,2'_Mmus_)[&&NHX:t=D])[&&NHX:t=S])[&&NHX:t=D],(3_Hsap_,(3_Mmus_,3'_Mmus_)[&&NHX:t=D])[&&NHX:t=S])[&&NHX:t=D];"
duplicates_tree_with_deletion = "(((1_Hsap_,1_Mmus_)[&&NHX:t=S],(2_Mmus_,2'_Mmus_)[&&NHX:t=D])[&&NHX:t=D],(3_Hsap_,(3_Mmus_,3'_Mmus_)[&&NHX:t=D])[&&NHX:t=S])[&&NHX:t=D];"
whiteboard_example_gene_tree = "(((a_Hsap_,a_Ptro_)[&&NHX:t=S],a_Mmus_)[&&NHX:t=S],(b_Hsap_,b_Mmus_)[&&NHX:t=S])[&&NHX:t=D];" 
whiteboard_example_species_tree = "((Hsap,Ptro),Mmus);"

st = Species_Tree(species_tree)
std = Species_Tree(species_tree)
wst = Species_Tree(whiteboard_example_species_tree)

gt = Gene_Tree(duplicates_tree)
gtd = Gene_Tree(duplicates_tree_with_deletion)
wgt = Gene_Tree(whiteboard_example_gene_tree)

st.add_gene_tree(gt)
st.print_count()


wst.print_count()
wst.add_gene_tree(wgt)
wst.print_count()


std.add_gene_tree(gtd)
std.print_count()

