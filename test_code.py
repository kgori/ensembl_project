#!/usr/bin/env python

import dendropy as dpy
from copy import copy,deepcopy
import itertools

species_tree = "(H,M);"
duplicates_tree = "(((H1,M1)[S],((H2,H2.)[D],(M2,M2.)[D])[S])[D],(H3,(M3,M3.)[D])[S])[D];"

st = dpy.Tree().get_from_string(species_tree,'newick')
dt = dpy.Tree().get_from_string(duplicates_tree,'newick',extract_comment_metadata=True)

# i = st.preorder_edge_iter()


i = list(st.preorder_edge_iter())
i[0].label=['alpha',0]
i[1].label=['beta',0]
i[2].label=['gamma',0]

index=0
s_edge = i[index]
# c=1 # just used for debugging
d = {}
if not s_edge.label: s_edge.label = 0
for d_edge in dt.preorder_edge_iter():
    # print "NODE ",c,"    ", # debug
    # c+=1 # debug
    try: 
        if 'D' in d_edge.tail_node.comments:
            if 'read' not in d_edge.tail_node.comments:
                s_edge.label[1] += 1
                d_edge.tail_node.comments.append('read')
                d[str(d_edge.tail_node)]=index
                # print "D found at unread tail node: {0}++ ({1})".format(s_edge.label[0], s_edge.label[1]) # debug
            else:
                index = d[str(d_edge.tail_node)]
                s_edge = i[index]
                # print "R! {0}= ({1})".format(s_edge.label[0], s_edge.label[1]) # debug


        elif 'S' in d_edge.tail_node.comments:
            try:
                index += 1
                s_edge = i[index]
                if not s_edge.label: s_edge.label = 1
                else: s_edge.label[1]+=1
                # print "S found at tail node: iter; {0}++ ({1})".format(s_edge.label[0], s_edge.label[1]) # debug
            except IndexError:
                print "Resetting...R! {0}= ({1})".format(s_edge.label[0], s_edge.label[1])
                i = st.preorder_edge_iter()
                s_edge = i.next()
                print s_edge.label

    except AttributeError: 
        s_edge.label[1] += 1
        # print "Root edge found: {0}++ ({1})".format(s_edge.label[0], s_edge.label[1]) # debug


for edge in st.preorder_edge_iter():
    print edge.label
