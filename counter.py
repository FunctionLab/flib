import logging
logger = logging.getLogger(__name__)

import sys
import numpy
import copy
import math
from cdatabase import CDatabase

class Counter:


    def __init__(self, cdata, quants=None):
        self.cdata = cdata

        if cdata is not None:
            self.datasets = cdata.get_datasets()

        self.global_std = None

        self.name = None
        self.prior = 0.01

        self.pos_counts = None
        self.neg_counts = None

        self.pos_pairs = None
        self.pos_pairs_cnt = None
        self.neg_pairs = None
        self.neg_pairs_cnt = None

        self.pos_cpt = None
        self.neg_cpt = None

        self.bin_effects = None

        self.quants = quants

    @classmethod
    def fromfilenames(cls, cdb_dir, didx, gidx, zeros_file, nibble=True, quants=None):
        cdata = CDatabase(cdb_dir, didx, gidx, zeros_file, nibble)
        return cls(cdata, quants)

    @classmethod
    def fromcountfile(cls, cfile, cdata=None):
        cfile = open(cfile)
        i = 0
        counter = cls(None, None)
        counter.cdata = cdata
        pos_counts = []
        neg_counts = []
        datasets = []
        for l in cfile:
            toks = l.strip().split('\t')
            if i == 0:
                (name, dcount) = toks
                counter.name = name
            elif i == 1:
                (gneg, gpos) = toks
                counter.pos_pairs_cnt = int(gpos)
                counter.neg_pairs_cnt = int(gneg)
                counter.prior = float(gpos)/(float(gneg)+float(gpos))
            elif (i+1) % 3 == 0:
                (dname,) = toks
                datasets.append(dname)
            elif (i+1) % 3 == 1:
                counts = [int(x) for x in toks]
                neg_counts.append(counts)
            elif (i+1) % 3 == 2:
                counts = [int(x) for x in toks]
                pos_counts.append(counts)
            i += 1

        if cdata is not None:
            counter.datasets = [None]*len(datasets)
            counter.pos_counts = [None]*len(datasets)
            counter.neg_counts = [None]*len(datasets)
            for (i,d) in enumerate(datasets):
                idx = cdata.get_dataset_idx(d)
                counter.datasets[idx] = d
                counter.pos_counts[idx] = pos_counts[i]
                counter.neg_counts[idx] = neg_counts[i]
        else:
            counter.datasets = datasets
            counter.pos_counts = pos_counts
            counter.neg_counts = neg_counts
            

        counter.neg_cpt = counter.calc_cpt(counter.neg_counts)
        counter.pos_cpt = counter.calc_cpt(counter.pos_counts)

        counter.bin_effects = counter.calc_bineffects()


        return counter

    def quantize(self, value):
        last = len(self.quants)-1
        if self.quants is not None:
            for i,q in enumerate(self.quants):
                if q >= value or i == last:
                    return i
        else:
            return value

    def set_global_std(self, std):
        self.global_std = std

    def get_cdata(self):
        return self.cdata

    def get_prior(self):
        return self.prior

    def get_counts(self):
        return (self.neg_counts, self.pos_counts)

    def get_bineffects(self):
        return self.bin_effects

    def get_cpt(self):
        return (self.neg_cpt, self.pos_cpt)

    def get_name(self):
        return self.name

    def get_datasets(self):
        return self.datasets[:]

    def get_bridge_negs(self, genes):
        if self.global_std is None:
            return None
        pairs = []
        glb_genes = self.global_std.get_genes()
        for g1 in genes:
            for (i,g2) in enumerate(glb_genes):

                if g2 in genes:
                    continue

                std = self.global_std.get_genepair_values(g1, g2)
                if len(std) > 0 and std[0] == 0:
                    pairs.append((g1, g2))

        return pairs


    def learn_ctxt(self, ctxt, genes, cpos=None, cneg=None, bpos=None, bneg=None):
        if (cpos is None) or (cneg is None) or (bpos is None) or (bneg is None):
            logger.error("Values must be set for learning: cpos:%s cneg:%s bpos:%s bneg:%s.", cpos, cneg, bpos, bneg)
            return False
        self.name = ctxt

        genes = set(genes)

        dbgenes = self.global_std.get_genes()

        MISSING = self.global_std.missing_val()
        pos_pairs = []
        neg_pairs = []
        for g in genes:
            std = self.global_std.get_gene_values(g)

            for (i, v) in enumerate(std):
                if v == MISSING:
                    continue
                other = dbgenes[i]
                if (other is None) or (other == g):
                    continue
                if other in genes:
                    if g > other:
                        continue
                    if v == 1:
                        if cpos:
                            pos_pairs.append((g, other))
                    else:
                        if cneg:
                            neg_pairs.append((g, other))
                else:
                    if v == 1:
                        if bpos:
                            pos_pairs.append((g, other))
                    else:
                        if bneg:
                            neg_pairs.append((g, other))

        self.pos_pairs = pos_pairs
        self.neg_pairs = neg_pairs 
        self.pos_pairs_cnt = len(pos_pairs)
        self.neg_pairs_cnt = len(neg_pairs)


        logger.info('Counting positives')
        self.pos_counts = self.count(self.pos_pairs)
        logger.warning('POS PAIRS: %s', len(self.pos_pairs))
        logger.info('Counting negatives')
        self.neg_counts = self.count(self.neg_pairs)
        logger.warning('NEG PAIRS: %s', len(self.neg_pairs))

        self.neg_cpt = self.calc_cpt(self.neg_counts)
        self.pos_cpt = self.calc_cpt(self.pos_counts)

        self.bin_effects = self.calc_bineffects()


    def count(self, gene_pairs):
        dsets = self.cdata.get_datasets()

        # Initialize with size of dsets (number of quants)
        dsets_counts = [None]*len(dsets)
        for i in range(len(dsets)):
            dsets_counts[i] = [0]*self.cdata.get_dataset_bins(i)

        MISSING = self.cdata.missing_val()
        for (i,pair) in enumerate(gene_pairs):

            v = self.cdata.get_genepair_values(pair[0], pair[1])
            for d,value in enumerate(v):
                if value == MISSING:
                    zbin = self.cdata.get_dataset_zero(d)
                    if zbin is not None:
                        value = zbin
                    else:
                        continue

                # Quantize if necessary
                if self.quants is not None:
                    value = self.quantize(value)

                dsets_counts[d][value] += 1

        return dsets_counts


    def calc_cpt(self,counts):
        cpt = copy.deepcopy(counts)
        #laplacian smoothing
        for i,d in enumerate(cpt):
            cpt[i] = [x+1 for x in d]

        for i in range(len(cpt)):
           total = float( sum(cpt[i]) ) # TODO: get the number of bins used
           for j in range(len(cpt[i])):
               if total > 0:
                   cpt[i][j] /= total
        return cpt

    def calc_bineffects(self):
        pos_log = []
        neg_log = []

        bin_effects = [None]*len(self.pos_cpt)
        for i in range(len(self.pos_cpt)):
            plog = numpy.log(self.pos_cpt[i])
            nlog = numpy.log(self.neg_cpt[i])
            bin_effects[i] = numpy.subtract(nlog,plog)

        return bin_effects


    def get_posterior(self, g1, g2, prior=None):
        v = self.cdata.get_genepair_values(g1, g2)

        if prior is None:
            prior = self.prior

        logratio = 0.0
        for i,value in enumerate(v):
            if self.cdata.is_missing(value):
                zbin = self.cdata.get_dataset_zero(i)
                if zbin is not None:
                    value = zbin
                else:
                    continue
            logratio += self.bin_effects[i][value]

        ratio = numpy.exp(logratio)*((1.0-float(prior))/float(prior))

        return 1.0/(1.0+ratio)


    def get_trust(self):
        datasets = self.datasets 
        for i in range(len(datasets)):
            trust = 0.0
            ratio = 0.0
            for j in range(len(self.pos_cpt[i])):
                pp = self.pos_cpt[i][j]
                pn = self.neg_cpt[i][j]
                diff = math.fabs(pp-pn)
                ratio += numpy.log(max(pp,pn) / min(pp,pn))
                trust += diff

            print self.name, datasets[i], trust, ratio
        return

    def print_counts(self):
        datasets = self.datasets
        print('\t'.join([self.name, str(len(datasets))]))
        print('\t'.join([str(self.neg_pairs_cnt), str(self.pos_pairs_cnt)]))
        for i,d in enumerate(datasets):
            print(d)
            print('\t'.join([str(x) for x in self.neg_counts[i]]))
            print('\t'.join([str(x) for x in self.pos_counts[i]]))


if __name__ == '__main__':
        from optparse import OptionParser

        usage = "usage: %prog [options]"
        parser = OptionParser(usage, version="%prog dev-unreleased")
        parser.add_option("-i", "--cdatabase-dir", dest="cdb", help="Directory of CDatabase", metavar="FILE")
        parser.add_option("-d", "--datasets", dest="dset", help="File of dataset names", metavar="FILE")
        parser.add_option("-g", "--gene-file", dest="gene_file", help="File of gene names", metavar="FILE")
        parser.add_option("-c", "--context-file", dest="context_file", help="File of gene names", metavar="FILE")
        parser.add_option("-z", "--zeros-file", dest="zeros_file", help="File of gene names", metavar="FILE")
        parser.add_option("-f", "--counts-file", dest="counts_file", help="Counts file", metavar="FILE")
        parser.add_option("-G", "--global-gold-standard", dest="gstd", help="global gold standard (CDatabase format)", metavar="FILE")
        parser.add_option("-l", "--gene1", dest="gene1", help="Query gene")
        parser.add_option("-L", "--gene2", dest="gene2", help="Query gene")
        parser.add_option("-q", "--ctxtpos", dest="cpos", help="Use positive edges between context genes", default=True)
        parser.add_option("-Q", "--ctxtneg", dest="cneg", help="Use negative edges between context genes", default=True)
        parser.add_option("-j", "--bridgepos", dest="bpos", help="Use bridging positives between context and non-context genes", default=False)
        parser.add_option("-J", "--bridgeneg", dest="bneg", help="Use bridging negatives between context and non-context genes", default=True)
        parser.add_option("-b", "--quant-file", dest="quant_file", help="File of bins for quantizing", metavar="FILE")
        parser.add_option("-B", "--byte", dest="byte", help="Size of data values", action="store_true", default=False)

        (options, args) = parser.parse_args()


        if options.cdb is None:
            sys.stderr.write("--cdatabase-dir is required.\n")
            sys.exit()
        if options.gene_file is None:
            sys.stderr.write("--gene-file is required.\n")
            sys.exit()

        genef = open(options.gene_file)
        gidx = []
        for l in genef:
            (idx, gene) = l.strip().split()
            gidx.append((int(idx)-1, gene))
        genef.close()

        dsf = open(options.dset)
        didx = []
        for l in dsf:
            (idx, ds, bins) = l.strip().split()
            didx.append((int(idx)-1, ds, int(bins)))
        dsf.close()

        quants = None
        if options.quant_file:
            quants = open(options.quant_file).readline().strip().split()
            quants = [float(x) for x in quants]


        # Load counter from a file or using a context file
        if options.counts_file:
            cdb = CDatabase(options.cdb, didx, gidx, options.zeros_file, not options.byte)
            counter = Counter.fromcountfile(options.counts_file, cdb)
        else:
            counter = Counter.fromfilenames(options.cdb, didx, gidx, options.zeros_file, not options.byte, quants)
            gold_std = None
            if options.gstd:
                gold_std = CDatabase(options.gstd, [(0, 'gold_standard', 2)], gidx)
                counter.set_global_std(gold_std)

            pos = None
            neg = None
            if options.context_file:
                genes_set = set()
                for l in open(options.context_file):
                    genes_set.add(l.strip())
                genes = list(genes_set)
                counts = counter.learn_ctxt(options.context_file, genes, cpos=options.cpos, cneg=options.cneg, bpos=options.bpos, bneg=options.bneg)




        qgenes = set()
        if options.gene1:
            qgenes.add(options.gene1)
        if options.gene2:
            qgenes.add(options.gene2)

        from network import Network
        net = Network.fromcounter(counter)

        print net.query(qgenes, .1)

