from optparse import OptionParser
import re
import sys
import os
from idmap import idmap
from go import go

# constants
# idx of reference IDs
ref_col = 5

usage = "usage: %prog [options]"
parser = OptionParser(usage, version = "%prog dev-unreleased")
parser.add_option("-o", "--obo-file", dest="obo", help="obo file", metavar="FILE")
parser.add_option("-a", "--association-file", dest="ass", help="gene association file", metavar="FILE")
parser.add_option("-b", dest="term_col", type="int", help = "What column of the annotations file contains the term identifiers?", default=4)
parser.add_option("-g", dest="gcol", type="int", help = "What column of the annotations file contains the desired identifiers?", default=1)
parser.add_option("-d", "--output-prefix", dest = "opref", help = "prefix for output files", metavar = "string")
parser.add_option("-f", "--output-filename", dest = "ofile", help = "If given outputs all go term/gene annotation pairs to this file, file is created in the output prefix directory.", metavar = "string")
parser.add_option("-i", "--id-file", dest = "idfile", help = "file to map excisting gene ids to the desired identifiers in the format <gene id>\\t<desired id>\\n", metavar = "FILE")
parser.add_option("-p", action="store_true", dest="progagate", help = "Should we progagate gene annotations?")
parser.add_option("-t", "--slim-file", dest="slim", help="GO slim file contains GO terms to output, if not given outputs all GO terms", metavar="FILE")
parser.add_option("-n", "--namespace", dest="nspace", help="limit the GO term output to the input namespace: (biological_process, cellular_component, molecular_function)", metavar="STRING")
parser.add_option("-r", dest="refids", action="store_true", help = "If given keeps track of ref IDs (e.g. PMIDs) for each go term and prints to standard out")
parser.add_option("-c", dest="check_fringe", action="store_true", help = "Is the given slim file a true fringe in the given obo file?")

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

if options.refids:
    gene_ontology.populate_annotations(options.ass, options.gcol, ref_col, term_col=options.term_col)
else:
    gene_ontology.populate_annotations(options.ass, options.gcol, term_col=options.term_col)

if options.idfile is not None:
    gene_ontology.map_genes(id_name)

if options.progagate:
    gene_ontology.propagate()

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
