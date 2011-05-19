import sys
import os
from idmap import idmap

class go:
    heads = None
    go_terms = None
    name2id = None
    id2gterm = None
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
        self.name2id = {}
        self.id2gterm = {}
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
                    gterm = GO_term(fields[1])
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

                self.name2id[gterm.name] = gterm.go_id
                self.id2gterm[gterm.go_id] = gterm.name
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
                    self.go_terms[pgo_id] = GO_term(pgo_id)

                gterm.is_a.append(self.go_terms[pgo_id])
                self.go_terms[pgo_id].parent_of.append(gterm)
            elif inside and fields[0] == 'relationship:':
                if fields[1].find('has_part') != -1:
                    #has part is not a parental relationship -- it is actually for children.
                    continue
                gterm.head = False
                pgo_id = fields[2]
                if not self.go_terms.has_key(pgo_id):
                    self.go_terms[pgo_id] = GO_term(pgo_id)
                gterm.relationship.append(self.go_terms[pgo_id])
                self.go_terms[pgo_id].parent_of.append(gterm)
            elif inside and fields[0] == 'is_obsolete:':
                gterm.head = False
                del self.go_terms[gterm.get_id()]
                del self.name2id[gterm.name]
                del self.id2gterm[gterm.go_id]


        #print "WTF!!", self.go_terms.has_key("GO:0010551")
    """
    propagate all gene annotations
    """
    def propagate(self):
        print >> sys.stderr, "Propagate gene annotations"
        for head_gterm in self.heads:
            print >> sys.stderr, head_gterm.name
            self.propagate_recurse(head_gterm)
    def propagate_recurse(self, gterm):
        if len(gterm.parent_of) == 0:
            return

        for child_term in gterm.parent_of:
            if child_term.namespace != gterm.namespace:
                continue
            self.propagate_recurse(child_term)

        # need to optimize
        for child_term in gterm.parent_of:
            if child_term.namespace != gterm.namespace:
                continue
            gterm.annotated_genes = gterm.annotated_genes | child_term.annotated_genes

            if len(child_term.refids.keys()) != 0:
                for refid in child_term.refids.keys():
                    if gterm.refids.has_key(refid):
                        gterm.refids[refid] = gterm.refids[refid] | child_term.refids[refid]
                    else:
                        gterm.refids[refid] = child_term.refids[refid].copy()


    def print_terms(self, out_dir, terms=None, p_namespace=None):
        if terms == None:
            terms = self.go_terms.keys()

        #print terms
        for gid in terms:
            if not self.go_terms.has_key(gid):
                if self.alt_id2std_id.has_key(gid):
                    gid = self.alt_id2std_id[gid]
                else:
                    print >> sys.stderr, "ERROR: term name don't exist:", gid
                    continue

            gterm = self.id2gterm[gid]

            if p_namespace != None and self.go_terms[gid].namespace != p_namespace:
                continue


            f = open(out_dir + '/' + gterm, 'w')

            for gene in self.go_terms[gid].annotated_genes:
                print >> f, gene
            f.close()

    def print_to_single_file(self, out_file, terms=None, p_namespace=None, gene_asso_format=False):
        f = open(out_file, 'w')
        if terms == None:
            terms = self.go_terms.keys()

        terms.sort()
        for gid in terms:
            if not self.go_terms.has_key(gid):
                if self.alt_id2std_id.has_key(gid):
                    gid = self.alt_id2std_id[gid]
                else:
                    print >> sys.stderr, "ERROR: term name don't exist:", gid
                    continue

            if p_namespace != None and self.go_terms[gid].namespace != p_namespace:
                continue

            for gene in self.go_terms[gid].annotated_genes:
                if gene_asso_format:
                    to_print = ['dummy', gene, 'dummy', 'dummy', gid]
                    print >> f, '\t'.join(to_print)
                else:
                    print >> f, gid + '\t' + gene
        f.close()

    # print each term ref IDs to a standard out
    def print_refids(self, terms=None, p_namespace=None):
        print >> sys.stderr, "printing ref IDs"

        if terms == None:
            terms = self.go_terms.keys()

        terms.sort()
        for gid in terms:
            if not self.go_terms.has_key(gid):
                if self.alt_id2std_id.has_key(gid):
                    gid = self.alt_id2std_id[gid]
                else:
                    print >> sys.stderr, "ERROR: term name don't exist:", gid
                    continue

            if p_namespace != None and self.go_terms[gid].namespace != p_namespace:
                continue

            for refid in self.go_terms[gid].refids.keys():
                for gene in self.go_terms[gid].refids[refid]:
                    print gid + '\t' + refid + '\t' + gene

    # be aware this is added only to be used with python sctript  cross_annotate_single_file_only_crossed.py
    def print_to_single_file_cross_annotated(self, out_file, terms=None, p_namespace=None):
        f = open(out_file, 'w')
        if terms == None:
            terms = self.go_terms.keys()

        terms.sort()
        for gid in terms:
            if not self.go_terms.has_key(gid):
                if self.alt_id2std_id.has_key(gid):
                    gid = self.alt_id2std_id[gid]
                else:
                    print >> sys.stderr, "ERROR: term name don't exist:", gid
                    continue

            if p_namespace != None and self.go_terms[gid].namespace != p_namespace:
                continue

            for gene in self.go_terms[gid].cross_annotated_genes:
                print >> f, gene + '\t' + gid
        f.close()

    def map_genes(self, id_name):
        for gid in self.go_terms.keys():
            self.go_terms[gid].map_genes(id_name)
            

    def populate_annotations(self, annotation_file, gene_col=1, ref_col=-1, term_col=4):
        f = open(annotation_file, 'r')
        print >> sys.stderr, "Populate gene annotations:", annotation_file
        for line in f:
            if line[0] == '!':
                continue
            fields = line.rstrip('\n').split('\t')
            go_term = fields[term_col]
            gene = fields[gene_col]

            if not self.go_terms.has_key(go_term):
                if self.alt_id2std_id.has_key(go_term):
                    go_term = self.alt_id2std_id[go_term]
                else:
                    print >> sys.stderr, "ERROR: term name don't exist:", go_term
                    continue

            self.go_terms[go_term].annotated_genes.add(gene)

            # if given a reference idx, store the reference ID (e.g. PMID)
            if ref_col != -1:
                refid = fields[ref_col]

                if self.go_terms[go_term].refids.has_key(refid):
                    self.go_terms[go_term].refids[refid].add(gene)
                else:
                    self.go_terms[go_term].refids[refid] = set([gene])

        f.close()
        self.populated = True


    def populate_additional_taxon_specificity(self, ncbi_tax_obj, taxon_specificity_add_file, tag_tax_id):
        f = open(taxon_specificity_add_file, 'r')
        print >> sys.stderr, "Populate GO specificity:", taxon_specificity_add_file

        if ncbi_tax_obj.id2species.has_key(tag_tax_id):
            self.go_organism_tax_id = tag_tax_id
        else:
            print >> sys.stderr, 'ERROR: NCBI tax ID do not exist', tag_tax_id
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
        f = open(taxon_specificity_obo_file, 'r')
        print >> sys.stderr, "Populate GO specificity:", taxon_specificity_obo_file

        if ncbi_tax_obj.id2species.has_key(tag_tax_id):
            self.go_organism_tax_id = tag_tax_id
        else:
            print >> sys.stderr, 'ERROR: NCBI tax ID do not exist', tag_tax_id
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
                            print >> sys.stderr, 'ERROR missing no NCBI tax id:', fl
                # now go label your go tree
                self.propagate_taxon_specificity(final_tax_ids, gid, relationship, ncbi_tax_obj)

                # ok now collected all info
                inside = False
                gid = None
                relationship = None
                tax_id = None

        f.close()

    def propagate_taxon_specificity(self, tax_ids, go_id, relationship, ncbi_tax_obj):
        current_gterm = None
        if self.go_terms.has_key(go_id):
            current_gterm = self.go_terms[go_id]
        elif self.alt_id2std_id.has_key(go_id):
            current_gterm = self.go_terms[self.alt_id2std_id[go_id]]
        else:
            print >> sys.stderr, "ERROR: term name don't exist:", go_id
            return


        if relationship == 'only_in_taxon':
            for tid in tax_ids:
                if ncbi_tax_obj.check_lineage(tid, self.go_organism_tax_id):
                    return
            self.propagate_taxon_set_false(go_id)
        elif relationship == 'never_in_taxon':
            for tid in tax_ids:
                if ncbi_tax_obj.check_lineage(tid, self.go_organism_tax_id):
                    self.propagate_taxon_set_false(go_id)
                    return
        else:
            print >> sys.stderr, "ERROR invalid relationship term", relationship
            return

    def propagate_taxon_set_false(self, gid):
        if not self.go_terms.has_key(gid):
            if  self.alt_id2std_id.has_key(gid):
                gid = self.alt_id2std_id[gid]
            else:
                print >> sys.stderr, "ERROR: term name don't exist:", gid
                return

        self.go_terms[gid].valid_go_term = False

        for child_term in self.go_terms[gid].parent_of:
            self.propagate_taxon_set_false(child_term.get_id())

    def goid2goterm(self, id):
        if self.id2gterm.has_key(id):
            return self.id2gterm[id]
        elif self.alt_id2std_id.has_key(id):
            return self.id2gterm[ self.alt_id2std_id[id] ]
        else:
            return ''
        
    # check if slim terms forms a true fringe in the obo structure
    def check_fringe(self, slim_file, namespace=None):
        leaf_goterms = []
        slim_goterms = []

        # add GO ids to the leaf terms
        for gid in self.go_terms.keys():
            if len(self.go_terms[gid].parent_of) == 0:
                if namespace != None and self.go_terms[gid].namespace != namespace:
                    continue

                self.go_terms[gid].annotated_genes.add(gid)
                leaf_goterms.append(gid)

        # now propagate the GO ids from the leaf terms
        self.propagate()

        # open go terms from slim term
        f = open(slim_file, 'r')
        gterms = []
        for line in f:
            fields = line.rstrip('\n').split('\t')
            gterms.append(fields[1])
        f.close()

        # print self.go_terms.keys()

        # now go colect the GO leaf term ids that have been propagated to the slim terms
        for gid in gterms:
            if not self.go_terms.has_key(gid):
                if self.alt_id2std_id.has_key(gid):
                    gid = self.alt_id2std_id[gid]
                else:
                    print >> sys.stderr, "ERROR: slim term name don't exist (probably obsolete term):", gid
                    continue

            slim_goterms.extend( self.go_terms[gid].annotated_genes )

        # now compare two sets
        leaf_goterms.sort()
        slim_goterms.sort()

        if leaf_goterms == slim_goterms:
            return True
        else:
            for lgoterm in leaf_goterms:
                if lgoterm not in slim_goterms:
                    print >> sys.stderr, "Missing leaf terms", self.id2gterm[lgoterm], lgoterm
            return False


class GO_term:
    go_id = ''
    is_a = None
    relationship = None
    parent_of = None
    annotated_genes = None
    alt_id = None
    namespace = ''
    included_in_all = None
    valid_go_term = None
    cross_annotated_genes = None
    head = None
    refids = None

    def __init__(self, go_id):
        self.head = True
        self.go_id = go_id
        self.annotated_genes = set([])
        self.refids = {}
        self.cross_annotated_genes = set([])
        self.is_a = []
        self.relationship = []
        self.parent_of = []
        self.alt_id = []
        self.included_in_all = True
        self.valid_go_term = True
    def get_id(self):
        return self.go_id
    def map_genes(self, id_name):
        mapped_genes_set = set([])
        for gene in self.annotated_genes:
            mapped_genes = id_name.get(gene)
            if mapped_genes == None:
                print >> sys.stderr, "ERROR: No matching gene id:", gene
                continue
            for mgene in mapped_genes:
                mapped_genes_set.add(mgene)
        self.annotated_genes = mapped_genes_set

        if len(self.refids.keys()) != 0:
            for refid in self.refids.keys():
                mapped_genes_set = set([])
                for gene in self.refids[refid]:
                    mapped_genes = id_name.get(gene)
                    if mapped_genes == None:
                        continue
                    for mgene in mapped_genes:
                        mapped_genes_set.add(mgene)
                self.refids[refid] = mapped_genes_set
