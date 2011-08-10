import sys
from counts import Counts
import logging
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)

import re
from idmap import idmap

anywhere = set()
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
                name = '_'.join(fields)
                name = re.sub('[^\w\s_-]', '_', name).strip().lower()
                name = re.sub('[-\s_]+', '_', name)
                gterm.name = name
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
                self.go_terms[pgo_id].parent_of.add(gterm)
                gterm.child_of.add(self.go_terms[pgo_id])
            elif inside and fields[0] == 'relationship:':
                if fields[1].find('has_part') != -1:
                    #has part is not a parental relationship -- it is actually for children.
                    continue
                gterm.head = False
                pgo_id = fields[2]
                if not self.go_terms.has_key(pgo_id):
                    self.go_terms[pgo_id] = GOTerm(pgo_id)
                gterm.relationship.append(self.go_terms[pgo_id])
                self.go_terms[pgo_id].parent_of.add(gterm)
                gterm.child_of.add(self.go_terms[pgo_id])
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
        if not len(gterm.parent_of):
            logger.debug("Base case with term %s", gterm.name)
            return

        for child_term in gterm.parent_of:
            self.propagate_recurse(child_term)
            new_annotations = set()
            for annotation in child_term.annotations:
                new_annotations.add(annotation.prop_copy())
            gterm.annotations = gterm.annotations | new_annotations

    """
    prune all gene annotations
    """
    def prune(self, eval_str):
        dterms = set()
        heads = set(self.heads)
        for (name, term) in self.go_terms.iteritems():
            if term in heads:
                print("Head term " + name)
                continue
            total = len(term.annotations)
            direct = 0
            for annotation in term.annotations:
                if annotation.direct:
                    direct += 1
            if eval(eval_str):
                for pterm in term.child_of:
                    pterm.parent_of.update(term.parent_of)
                    pterm.parent_of.discard(term)
                for cterm in term.parent_of:
                    cterm.child_of.update(term.child_of)
                    cterm.child_of.discard(term)
                dterms.add(name)
        for name in dterms:
            del self.go_terms[name]

        #remove connections to root if there are other parents
        for (name, term) in self.go_terms.iteritems():
            #if there is something in the intersection
            intersection = term.child_of & heads
            if (intersection):
                #if the intersection isn't the only thing it's a child of
                if (term.child_of - heads):
                    term.child_of -= intersection
                    for hterm in intersection:
                        hterm.parent_of.remove(term)

    def get_term(self, tid):
        logger.debug('get_term: %s', tid)
        term = None
        try:
            term = self.go_terms[tid]
        except KeyError:
            try:
                term = self.go_terms[self.alt_id2std_id[tid]]
            except KeyError:
                logger.error('Term name does not exist: %s', tid)
        return term

    def get_termobject_list(self, terms=None, p_namespace=None):
        logger.info('get_termobject_list')
        if terms is None:
            terms = self.go_terms.keys()
        reterms = []
        for tid in terms:
            obo_term = self.get_term(tid)
            if obo_term is None:
                continue
            if p_namespace is not None and obo_term.namespace != p_namespace:
                continue
            reterms.append(obo_term)
        return reterms

    def get_termdict_list(self, terms=None, p_namespace=None):
        logger.info('get_termdict_list')
        tlist = self.get_termobject_list(terms=terms, p_namespace=p_namespace)
        reterms = []
        for obo_term in tlist:
            reterms.append({'oboid':obo_term.go_id, 'name':obo_term.name})
        return reterms

    def print_terms(self, out_dir, terms=None, p_namespace=None):
        logger.info('print_terms')
        tlist = self.get_termobject_list(terms=terms, p_namespace=p_namespace)
        #print terms
        for term in tlist:
            f = open(out_dir + '/' + term.name, 'w')
            for annotation in term.annotations:
                print >> f, annotation.gid
            f.close()

    def print_to_single_file(self, out_file, terms=None, p_namespace=None, gene_asso_format=False):
        logger.info('print_to_single_file')
        tlist = self.get_termobject_list(terms=terms, p_namespace=p_namespace)
        tlist.sort()
        f = open(out_file, 'w')
        for term in tlist:
            for annotation in term.annotations:
                if gene_asso_format:
                    to_print = [annotation.xdb if annotation.xdb else '',
                                annotation.gid if annotation.xdb else '',
                                '', '', #Gene Symbol, NOT/''
                                term.go_id if term.go_id else '',
                                annotation.ref if annotation.ref else '',
                                annotation.evidence if annotation.evidence else '',
                                annotation.date if annotation.date else '',
                                str(annotation.direct)] #Direct is added in to indicate prop status
                    print >> f, '\t'.join([str(x) for x in to_print])
                else:
                    print >> f, term.go_id + '\t' + annotation.gid
        f.close()


    def dictify(self, term, thedict):
        direct = 0
        total = len(term.annotations)
        for annotation in term.annotations:
            if annotation.direct:
                direct += 1
        child_vals = []
        for child in term.parent_of:
            cdict = {}
            self.dictify(child, cdict)
            child_vals.append(cdict)
        thedict["name"] = term.name
        thedict["direct"] = direct
        thedict["total"] = total
        if child_vals:
            thedict["children"] = child_vals
        return

    def to_json(self):
        """
        Return the hierarchy for all nodes with more than min genes
        as a json string (depends on simplejson).
        """
        import simplejson
        redict = {}
        for head in self.heads:
            self.dictify(head, redict)
        return 'var ontology = ' + simplejson.dumps(redict, indent=2)

    # print each term ref IDs to a standard out
    def print_refids(self, terms=None, p_namespace=None):
        logger.info('print_refids')
        tlist = self.get_termobject_list(terms=terms, p_namespace=p_namespace)
        tlist.sort()
        for term in tlist:
            for annotation in term.annotations:
                print term.go_id + '\t' + annotation.ref + '\t' + annotation.gid

    # be aware this is added only to be used with python script  cross_annotate_single_file_only_crossed.py
    def print_to_single_file_cross_annotated(self, out_file, terms=None, p_namespace=None):
        logger.info('print_to_single_file_cross_annotated')
        tlist = self.get_termobject_list(terms=terms, p_namespace=p_namespace)
        tlist.sort()
        f = open(out_file, 'w')
        for term in tlist:
            for gene in term.cross_annotated_genes:
                print >> f, gene + '\t' + term.go_id
        f.close()

    def map_genes(self, id_name):
        for go_term in self.go_terms.itervalues():
            go_term.map_genes(id_name)

    def populate_annotations(self, annotation_file, xdb_col=0, gene_col=None, term_col=None, ref_col=5, ev_col=6, date_col=13):
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
            try:
                ref = fields[ref_col]
            except IndexError:
                ref = None
            try:
                ev = fields[ev_col]
            except IndexError:
                ev = None
            try:
                date = fields[date_col]
            except IndexError:
                date = None

            try:
                details = fields[details_col]
                if details == 'NOT':
                    continue
            except IndexError:
                pass
            go_term = self.get_term(go_id)
            if go_term is None:
                continue
            logger.info('Gene %s and term %s', gene, go_term.go_id)
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
    child_of = None
    annotations = None
    alt_id = None
    namespace = ''
    included_in_all = None
    valid_go_term = None
    cross_annotated_genes = None
    head = None
    name = None
    base_counts = None
    counts = None
    weight = None

    def __init__(self, go_id):
        self.head = True
        self.go_id = go_id
        self.annotations = set([])
        self.cross_annotated_genes = set([])
        self.is_a = []
        self.relationship = []
        self.parent_of = set()
        self.child_of = set()
        self.alt_id = []
        self.included_in_all = True
        self.valid_go_term = True
        self.name = None
        self.base_counts = Counts()
        self.counts = Counts()
        self.weight = 0.0

    def __cmp__(self, other):
        return cmp(self.go_id, other.go_id)

    def __hash__(self):
        return(self.name.__hash__())

    def read_base_counts(self, filename):
        self.base_counts.read(filename)

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

    def get_annotated_genes(self):
        genes = []
        for annotation in self.annotations:
            genes.append(annotation.gid)
        return genes

    def add_annotation(self, gid):
        self.annotations.add(Annotation(gid=gid))

if __name__ == '__main__':
    from optparse import OptionParser

    usage = "usage: %prog [options]"
    parser = OptionParser(usage, version="%prog dev-unreleased")
    parser.add_option("-o", "--obo-file", dest="obo", help="obo file", metavar="FILE")
    parser.add_option("-a", "--association-file", dest="ass", help="gene association file", metavar="FILE")
    parser.add_option("-b", dest="term_col", type="int", help="What column of the annotations file contains the term identifiers?", default=4)
    parser.add_option("-g", dest="gcol", type="int", help="What column of the annotations file contains the desired identifiers?", default=1)
    parser.add_option("-d", "--output-prefix", dest="opref", help="prefix for output files", metavar="string")
    parser.add_option("-f", "--output-filename", dest="ofile", help="If given outputs all go term/gene annotation pairs to this file, file is created in the output prefix directory.", metavar="string")
    parser.add_option("-i", "--id-file", dest="idfile", help="file to map excisting gene ids to the desired identifiers in the format <gene id>\\t<desired id>\\n", metavar="FILE")
    parser.add_option("-p", action="store_true", dest="progagate", help="Should we progagate gene annotations?")
    parser.add_option("-P", "--prune", dest="prune", help="A python string that will be evaled to decide if a node should be pruned.  Available variables are 'total' and 'direct' which are the total number of annotations and the number of direct annotations.")
    parser.add_option("-t", "--slim-file", dest="slim", help="GO slim file contains GO terms to output, if not given outputs all GO terms", metavar="FILE")
    parser.add_option("-n", "--namespace", dest="nspace", help="limit the GO term output to the input namespace: (biological_process, cellular_component, molecular_function)", metavar="STRING")
    parser.add_option("-r", dest="refids", action="store_true", help="If given keeps track of ref IDs (e.g. PMIDs) for each go term and prints to standard out")
    parser.add_option("-c", dest="check_fringe", action="store_true", help="Is the given slim file a true fringe in the given obo file?  Prints the result and exits.")
    parser.add_option("-j", "--json-file", dest="json", help="file to output ontology (as json) to.")
    (options, args) = parser.parse_args()

    if options.obo is None:
        sys.stderr.write("--obo file is required.\n")
        sys.exit()
    if options.check_fringe is None and options.ass is None:
        sys.stderr.write("--association file is required.\n")
        sys.exit()
    if options.check_fringe is None and options.opref is None and not options.refids:
        sys.stderr.write("--prefix is required.\n")
        sys.exit()
    if options.check_fringe and options.slim is None:
        sys.stderr.write("--When checking fringe, must provide slim file.\n")
        sys.exit()

    id_name = None
    if options.idfile is not None:
        id_name = idmap(options.idfile)

    gene_ontology = go(options.obo)

    # only check if fringe is valid in this obo file?
    if options.check_fringe:
        if gene_ontology.check_fringe(options.slim, options.nspace):
            print "A complete fringe"
        else:
            print "not a fringe"
        # now exit
        sys.exit(0)

    gene_ontology.populate_annotations(options.ass, gene_col=options.gcol, term_col=options.term_col)

    if options.idfile is not None:
        gene_ontology.map_genes(id_name)

    if options.progagate:
        gene_ontology.propagate()

    if options.prune:
        gene_ontology.prune(options.prune)

    if options.json:
        jsonstr = gene_ontology.to_json()
        f = open(options.json, 'w')
        f.write(jsonstr)
        f.close()

    if options.slim:
        f = open(options.slim, 'r')
        gterms = []
        for line in f:
            fields = line.rstrip('\n').split('\t')
            gterms.append(fields[1])
        f.close()

        # should I only print ref IDs?
        if options.refids:
            gene_ontology.print_refids(gterms, options.nspace)
        elif options.ofile:
            gene_ontology.print_to_single_file(options.opref + '/' + options.ofile, gterms, options.nspace)
        else:
            gene_ontology.print_terms(options.opref, gterms, options.nspace)
    else:
        if options.refids:
            gene_ontology.print_refids(None, options.nspace)
        elif options.ofile:
            gene_ontology.print_to_single_file(options.opref + '/' + options.ofile, None, options.nspace)
        else:
            gene_ontology.print_terms(options.opref, None, options.nspace)

