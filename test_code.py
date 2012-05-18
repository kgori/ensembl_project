#!/usr/bin/python
# -*- coding: utf-8 -*-

from genecounter import GeneTree, SpeciesTree, TreeAlignmentError, \
    NodeLabellingError

import sys

if len(sys.argv) > 1:
    st = SpeciesTree(open(sys.argv[1]).read(), format=1)
    gt = GeneTree(open(sys.argv[2]).read(), format=1)
    if not '[&&NHX' in open(sys.argv[2]).read():
        gt.get_descendant_evol_events()
else:

    st = \
        SpeciesTree('(SE001:60.971176,(SE002:54.097175,((SE003:8.863838,SE013:8.7271048):33.859374,(((SE004:40.303116,(SE007:35.957807,(SE008:27.55188,(SE010:23.561527,SE011:23.148886):2.2511886):0.10991615):3.753085):7.938699,SE006:35.196556):7.992514,((SE005:16.104527,(SE014:1.9441679,SE015:3.2796744):19.411262):15.863363,(SE009:19.907848,SE012:39.749202):5.4071946):19.687332):3.7807204):17.137093):7.2302568):0;'
                    )
    gt = \
        GeneTree('(SE001,(SE002,((((SE004,(SE007,(SE008,SE010))),SE006),((SE004,(SE007,(SE008,(SE008,SE010))[&&NHX:D=Y])),SE006))[&&NHX:D=Y],(((SE005,(SE014,SE015)),(SE009,SE012)),(SE005,(SE009,SE012)))[&&NHX:D=Y])));'
                 , format=1)
for node in gt.traverse():
    node.identify()
for node in st.traverse():
    node.identify()

# gt.display()
gt._get_species_subtree(st, verbose=True)
st._finalise_counts()
for node in gt.traverse():
    if node.ambiguous:
        print "AMBIG",node
st.display_with_gene_tree(gt)
