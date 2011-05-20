import sys
import logging
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

import os
from idmap import idmap

class go:
    heads = None
    go_terms = None
    alt_id2std_id = None
    populated = None

    # populate this field if you want to mark this GO as organism specific
    go_organism_tax_id = None

    """
    Pass the obo file
    """
    def __init__(self, obo_file):
        self.heads = []
        self.go_terms = {}
        self.alt_id2std_id = {}
        self.populated = False

        f = open(obo_file, 'r')
        inside = False
        gterm = None
        for line in f:
            fields = line.rstrip().split()

            if len(fields) < 1:
                continue
            elif fields[0] == '[Term]':
                if gterm:
                    if gterm.head:
                        self.heads.append(gterm)
                inside = True
            elif fields[0] == '[Typedef]':
                if gterm:
                    if gterm.head:
                        self.heads.append(gterm)
                inside = False

            elif inside and fields[0] == 'id:':
                #print fields[1]
                if self.go_terms.has_key(fields[1]):
                    gterm = self.go_terms[fields[1]]
                else:
                    gterm = GOTerm(fields[1])
                    self.go_terms[ gterm.get_id() ] = gterm
                #print self.go_terms[fields[1]]
            elif inside and fields[0] == 'name:':
                fields.pop(0)
                gterm.name = '_'.join(fields)
                gterm.name = gterm.name.replace("'", "_")
                gterm.name = gterm.name.replace("-", "_")
                gterm.name = gterm.name.replace(",", "_")
                gterm.name = gterm.name.replace("/", "_")
                gterm.name = gterm.name.replace("+", "_")
                gterm.name = gterm.name.replace("(", "_")
                gterm.name = gterm.name.replace(")", "_")
            elif inside and fields[0] == 'namespace:':
                gterm.namespace = fields[1]
            elif inside and fields[0] == 'alt_id:':
                gterm.alt_id.append(fields[1])
                self.alt_id2std_id[fields[1]] = gterm.get_id()
            elif inside and fields[0] == 'is_a:':
                gterm.head = False
                fields.pop(0)
                pgo_id = fields.pop(0)
                if not self.go_terms.has_key(pgo_id):
                    self.go_terms[pgo_id] = GOTerm(pgo_id)

                gterm.is_a.append(self.go_terms[pgo_id])
                self.go_terms[pgo_id].parent_of.append(gterm)
            elif inside and fields[0] == 'relationship:':
                if fields[1].find('has_part') != -1:
                    #has part is not a parental relationship -- it is actually for children.
                    continue
                gterm.head = False
                pgo_id = fields[2]
                if not self.go_terms.has_key(pgo_id):
                    self.go_terms[pgo_id] = GOTerm(pgo_id)
                gterm.relationship.append(self.go_terms[pgo_id])
                self.go_terms[pgo_id].parent_of.append(gterm)
            elif inside and fields[0] == 'is_obsolete:':
                gterm.head = False
                del self.go_terms[gterm.get_id()]

    """
    propagate all gene annotations
    """
    def propagate(self):
        logger.info("Propagate gene annotations")
        for head_gterm in self.heads:
            logger.info("Propagating %s", head_gterm.name)
            self.propagate_recurse(head_gterm)

    def propagate_recurse(self, gterm):
        if len(gterm.parent_of) == 0:
            logger.debug("Base case with term %s", gterm.name)
            return

        for child_term in gterm.parent_of:
            if child_term.namespace != gterm.namespace:
                continue
            self.propagate_recurse(child_term)
            new_annotations = set()
            for annotation in child_term.annotations:
                new_annotations.add(annotation.prop_copy())
            gterm.annotations = gterm.annotations | new_annotations

    def get_term(self, tid):
        term = None
        try:
            term = self.go_terms[tid]
        except KeyError:
            try:
                term = self.go_terms[self.alt_id2std_id[tid]]
            except KeyError:
                logger.error('Term name does not exist: %s', tid)
        return term

    def print_terms(self, out_dir, terms=None, p_namespace=None):
        logger.info('Print terms')
        if terms == None:
            terms = self.go_terms.keys()

        #print terms
        for tid in terms:
            go_term = self.get_term(tid)
            if go_term is None:
                continue

            if p_namespace != None and go_term.namespace != p_namespace:
                continue

            f = open(out_dir + '/' + go_term.name, 'w')
            for annotation in go_term.annotations:
                print >> f, annotation.gid
            f.close()

    def print_to_single_file(self, out_file, terms=None, p_namespace=None, gene_asso_format=False):
        logger.info('Printing to single file')
        f = open(out_file, 'w')
        if terms == None:
            terms = self.go_terms.keys()

        terms.sort()
        for tid in terms:
            go_term = self.get_term(tid)
            if go_term is None:
                continue
            if p_namespace != None and go_term.namespace != p_namespace:
                continue

            for annotation in go_term.annotations:
                if gene_asso_format:
                    to_print = [annotation.xdb if annotation.xdb else '',
                                annotation.gid,
                                '', '', #Gene Symbol, NOT/''
                                tid,
                                annotation.ref,
                                annotation.evidence,
                                annotation.date,
                                annotation.direct] #Direct is added in to indicate prop status
                    print >> f, '\t'.join(to_print)
                else:
                    print >> f, tid + '\t' + annotation.gid
        f.close()

    # print each term ref IDs to a standard out
    def print_refids(self, terms=None, p_namespace=None):
        logger.info('Printing ref IDs')

        if terms == None:
            terms = self.go_terms.keys()

        terms.sort()
        for tid in terms:
            go_term = self.get_term(gid)
            if go_term is None:
                continue
            if p_namespace != None and go_term.namespace != p_namespace:
                continue

            for annotation in go_term.annotations:
                print tid + '\t' + annotation.ref + '\t' + annotation.gid

    # be aware this is added only to be used with python script  cross_annotate_single_file_only_crossed.py
    def print_to_single_file_cross_annotated(self, out_file, terms=None, p_namespace=None):
        logger.info('Printing to single file, cross annotated')
        f = open(out_file, 'w')
        if terms == None:
            terms = self.go_terms.keys()

        terms.sort()
        for tid in terms:
            go_term = self.get_term(tid)
            if go_term is None:
                continue

            if p_namespace != None and go_term.namespace != p_namespace:
                continue

            for gene in go_term.cross_annotated_genes:
                print >> f, gene + '\t' + tid
        f.close()

    def map_genes(self, id_name):
        for go_term in self.go_terms.itervalues():
            go_term.map_genes(id_name)


    def populate_annotations(self, annotation_file, xdb_col=0, gene_col=1, term_col=4, ref_col=5, ev_col=6, date_col=13):
        logger.info('Populate gene annotations: %s', annotation_file)
        details_col = 3
        f = open(annotation_file, 'r')
        for line in f:
            if line[0] == '!':
                continue
            fields = line.rstrip('\n').split('\t')
            xdb = fields[xdb_col]
            gene = fields[gene_col]
            go_id = fields[term_col]
            ref = fields[ref_col]
            ev = fields[ev_col]
            date = fields[date_col]

            details = fields[details_col]
            if details == 'NOT':
                continue

            go_term = self.get_term(go_id)
            if go_term is None:
                continue
            annotation = Annotation(xdb=xdb, gid=gene, ref=ref, evidence=ev, date=date, direct=True)
            go_term.annotations.add(annotation)

        f.close()
        self.populated = True


    def populate_additional_taxon_specificity(self, ncbi_tax_obj, taxon_specificity_add_file, tag_tax_id):
        logger.info("Populate GO specificity: %s", taxon_specificity_add_file)

        f = open(taxon_specificity_add_file, 'r')

        if ncbi_tax_obj.id2species.has_key(tag_tax_id):
            self.go_organism_tax_id = tag_tax_id
        else:
            logger.error("NCBI tax ID %s does not exist", tag_tax_id)
            sys.exit(1)

        for line in f:
            fields = line.rstrip('\n').split('\t')
            if len(fields) == 0:
                continue

            if line[0] == '#':
                continue

            gid = fields[0]
            relationship = fields[1]
            org = fields[2]
            # now go label your go tree
            self.propagate_taxon_specificity([org], gid, relationship, ncbi_tax_obj)

        f.close()


    def populate_taxon_specificity(self, ncbi_tax_obj, taxon_specificity_obo_file, tag_tax_id):
        logger.info("Populate GO specificity: %s", taxon_specificity_obo_file)
        f = open(taxon_specificity_obo_file, 'r')
        if ncbi_tax_obj.id2species.has_key(tag_tax_id):
            self.go_organism_tax_id = tag_tax_id
        else:
            logger.error("NCBI tax ID %s does not exist", tag_tax_id)
            sys.exit(1)

        inside = False
        gid = None
        relationship = None
        tax_id = None

        only_in_taxon = set([])
        never_in_taxon = set([])
        for line in f:
            fields = line.rstrip('\n').split()
            if len(fields) == 0:
                continue

            if fields[0] == '[Term]':
                inside = True
            elif inside and fields[0] == 'id:':
                gid = fields[1]
            elif inside and fields[0] == 'relationship:':
                relationship = fields[1]
                (tax_type, tax_id) = fields[2].split(':')

                final_tax_ids = []
                if tax_type == 'NCBITaxon':
                    final_tax_ids.append(tax_id)
                elif tax_type == 'NCBITaxon_Union':
                    for i, fl in enumerate(fields):
                        if not (i > 3 and fl != 'or'):
                            continue

                        if ncbi_tax_obj.species2id.has_key(fl):
                            final_tax_ids.append(ncbi_tax_obj.species2id[fl])
                        elif ncbi_tax_obj.in_part.has_key(fl):
                            [ final_tax_ids.append(in_part_id) for in_part_id in ncbi_tax_obj.in_part[fl] ]
                        else:
                            logger.error("Missing NCBI tax ID: %s", fl)
                # now go label your go tree
                self.propagate_taxon_specificity(final_tax_ids, gid, relationship, ncbi_tax_obj)

                # ok now collected all info
                inside = False
                gid = None
                relationship = None
                tax_id = None

        f.close()

    def propagate_taxon_specificity(self, tax_ids, term_id, relationship, ncbi_tax_obj):
        current_gterm = self.get_term(term_id)
        if current_gterm is None:
            return


        if relationship == 'only_in_taxon':
            for tid in tax_ids:
                if ncbi_tax_obj.check_lineage(tid, self.go_organism_tax_id):
                    return
            self.propagate_taxon_set_false(term_id)
        elif relationship == 'never_in_taxon':
            for tid in tax_ids:
                if ncbi_tax_obj.check_lineage(tid, self.go_organism_tax_id):
                    self.propagate_taxon_set_false(term_id)
                    return
        else:
            logger.error('Invalid relationship term: %s', relationship)
            return

    def propagate_taxon_set_false(self, tid):
        go_term = self.get_term(tid)
        if go_term is None:
            return

        go_term.valid_go_term = False

        for child_term in go_term.parent_of:
            self.propagate_taxon_set_false(child_term.get_id())

    # check if slim terms forms a true fringe in the obo structure
    def check_fringe(self, slim_file, namespace=None):
        leaf_tids = []
        slim_tids = []

        # add GO ids to the leaf terms
        for tid in self.go_terms.keys():
            leaf_term = self.go_terms[tid]
            if len(leaf_term.parent_of) == 0:
                if namespace != None and leaf_term.namespace != namespace:
                    continue
                leaf_term.annotations.add(Annotation(gid=tid))
                leaf_tids.append(tid)

        # now propagate the GO ids from the leaf terms
        self.propagate()

        # open go terms from slim term
        f = open(slim_file, 'r')
        stids = []
        for line in f:
            fields = line.rstrip('\n').split('\t')
            stids.append(fields[1])
        f.close()

        # now go colect the GO leaf term ids that have been propagated to the slim terms
        for tid in stids:
            slim_term = self.get_term(tid)
            if slim_term is None:
                logger.error('Slim term name does not exist (potentially obsolete term): %s', gid)
                continue
            slim_tids.extend([annotation.gid for annotation in slim_term.annotations])

        # now compare two sets
        leaf_tids.sort()
        slim_tids.sort()

        if leaf_tids == slim_tids:
            return True
        else:
            for lgoterm in leaf_tids:
                if lgoterm not in slim_tids:
                    logger.warning("Missing leaf terms: %s", lgoterm)
            return False

class Annotation(object):
    def __init__(self, xdb=None, gid=None, ref=None, evidence=None, date=None, direct=False):
        super(Annotation, self).__setattr__('xdb', xdb)
        super(Annotation, self).__setattr__('gid', gid)
        super(Annotation, self).__setattr__('ref', ref)
        super(Annotation, self).__setattr__('evidence', evidence)
        super(Annotation, self).__setattr__('date', date)
        super(Annotation, self).__setattr__('direct', direct)

    def prop_copy(self):
        return Annotation(xdb=self.xdb, gid=self.gid, ref=self.ref,
                          evidence=self.evidence, date=self.date, direct=False)

    def __hash__(self):
        return hash((self.xdb, self.gid, self.ref, self.evidence,
                     self.date, self.direct))

    def __eq__(self, other):
        return (self.xdb, self.gid, self.ref, self.evidence, self.date,
                self.direct).__eq__((other.xdb, other.gid, other.ref,
                                     other.evidence, other.date, other.direct))

    def __setattr__(self, *args):
        raise TypeError("Attempt to modify immutable object.")
    __delattr__ = __setattr__

class GOTerm:
    go_id = ''
    is_a = None
    relationship = None
    parent_of = None
    annotations = None
    alt_id = None
    namespace = ''
    included_in_all = None
    valid_go_term = None
    cross_annotated_genes = None
    head = None
    name = None

    def __init__(self, go_id):
        self.head = True
        self.go_id = go_id
        self.annotations = set([])
        self.cross_annotated_genes = set([])
        self.is_a = []
        self.relationship = []
        self.parent_of = []
        self.alt_id = []
        self.included_in_all = True
        self.valid_go_term = True
        self.name = None
    def get_id(self):
        return self.go_id
    def map_genes(self, id_name):
        mapped_annotations_set = set([])
        for annotation in self.annotations:
            mapped_genes = id_name.get(annotation.gid)
            if mapped_genes == None:
                logger.warning('No matching gene id: %s', annotation.gid)
                continue
            for mgene in mapped_genes:
                mapped_annotations_set.add(Annotation(xdb=None, gid=mgene,
                                                direct=annotation.direct,
                                                ref=annotation.ref,
                                                evidence=annotation.evidence,
                                                date=annotation.date))
        self.annotations = mapped_annotations_set

