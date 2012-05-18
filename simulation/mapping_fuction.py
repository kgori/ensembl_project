#!/usr/bin/python
# -*- coding: utf-8 -*-


def _get_species_subtree(self, species_tree):

    # define terms

    gtp = self  # gtp = GeneTree pointer, i.e. pointer to a node on the gene tree
    stp = species_tree  # stp = SpeciesTree pointer
    gtp.identify()
    stp.identify()

    if not stp.leaf:
        stp.L = species_tree.children[0]
        stp.R = species_tree.children[1]

    if gtp.leaf:
        if stp.leaf:
            assert gtp._taxon_set == stp._taxon_set
            gtp._species_subtree = stp
            if gtp.root:
                stp._count += 1
            return
        elif gtp._taxon_set <= stp.L._taxon_set:
            if not gtp.root:
                stp.R._losses += 1
            gtp._get_species_subtree(stp.L)
            return
        elif gtp._taxon_set <= stp.R._taxon_set:
            if not gtp.root:
                stp.L._losses += 1
            gtp._get_species_subtree(stp.R, verbose)
            return
        else:
            raise TreeAlignmentError('Gene tree leaf doesn\'t match species tree'
                    )
    else:
        gtp.L = self.children[0]
        gtp.R = self.children[1]

        if stp.leaf:
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
            gtp._species_subtree = stp
            stp._gains += 1
            if gtp.root:
                stp._count += 1
            gtp.L._get_species_subtree(stp, verbose)
            gtp.R._get_species_subtree(stp, verbose)
            return
        elif gtp.L._taxon_set <= stp.L._taxon_set and gtp.L._taxon_set \
            <= stp.R._taxon_set:
            raise TreeAlignmentError('Possible Mislabelling of GeneTree nodes'
                    )
        elif gtp.L._taxon_set <= stp.L._taxon_set and gtp.R._taxon_set \
            <= stp.L._taxon_set:
            if not gtp.root:
                stp.R._losses += 1
            gtp._get_species_subtree(stp.L, verbose)
            return
        elif gtp.L._taxon_set <= stp.L._taxon_set and gtp.R._taxon_set \
            <= stp.R._taxon_set:
            gtp._species_subtree = stp
            if gtp.root:
                stp._count += 1
            if gtp.duplication:
                raise TreeAlignmentError('shouldn\'t be a duplication here...'
                        )
            gtp.L._get_species_subtree(stp.L, verbose)
            gtp.R._get_species_subtree(stp.R, verbose)
            return
        elif gtp.L._taxon_set <= stp.R._taxon_set and gtp.R._taxon_set \
            <= stp.L._taxon_set:
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
        elif gtp.L._taxon_set <= stp.R._taxon_set and gtp.R._taxon_set \
            <= stp.R._taxon_set:
            if not gtp.root:
                stp.L._losses += 1
            gtp._get_species_subtree(stp.R, verbose)
            return
        elif gtp.R._taxon_set <= stp.R._taxon_set and gtp.R._taxon_set \
            <= stp.L._taxon_set:
            raise TreeAlignmentError('Possible Mislabelling of GeneTree nodes'
                    )
        elif stp.L._taxon_set & gtp._taxon_set and stp.R._taxon_set \
            & gtp._taxon_set:
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

            stp._gains += 1
            if gtp.root:
                stp._count += 1
            gtp._species_subtree = stp
            gtp.L._get_species_subtree(stp, verbose)
            gtp.R._get_species_subtree(stp, verbose)
            return
        else:
            raise TreeAlignmentError('Case fell through')
    return
