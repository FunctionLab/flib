import logging
logger = logging.getLogger(__name__)

import operator
import sys
import numpy
from counter import Counter

class Network:

    def __init__(self, bineffects, cdata, prior):
        self.bineffects = bineffects
        self.cdata = cdata
        self.genes = cdata.get_genes()
        self.datasets = cdata.get_datasets()
        self.prior = prior
        self.MISSING = cdata.missing_val()
        self.logprior = numpy.log((1.0-self.prior)/self.prior)

    @classmethod
    def fromcounter(cls, counter):
        return cls(counter.get_bineffects(), counter.get_cdata(), counter.get_prior())


    """
    Convenience method to calculate posterior, given a list of bin values
    """
    def calc_posterior(self, bin_values):
        logratio = 0.0
        for (i,b) in enumerate(bin_values):
            if b == self.MISSING:
                continue
            logratio += self.bineffects[i][b]
        logratio += self.logprior
        return self.to_posterior(logratio)

    """
    Convenience method to calculate posterior of a gene pair
    """
    def get_posterior(self, g1, g2):
        v = self.cdata.get_genepair_values(g1, g2)
        return self.calc_posterior(v)

    """
    Convenience method converting PR(d|^FR)/PR(d|FR) to PR(FR|D)
    """
    def to_posterior(self, logratio):
        ratio = numpy.exp(logratio)
        return 1.0/(1.0+ratio)

    """
    Return a list of edges (as tuples) between all gene pairs in given gene set
    """
    def get_edges(self, genes, edge_weight = 0, cache = None):
        genes_list = list(genes)
        edges = []
        for i in range(len(genes_list)):
            g1 = genes_list[i]
            for j in range(i+1,len(genes_list)):
                g2 = genes_list[j]
                posterior = self.get_posterior(g1, g2)
                if posterior >  edge_weight:
                    edges.append((g1, g2, posterior))
        return edges

    """
    Return the top datasets contributing to an edge posterior
    """
    def get_evidence(self, g1, g2, top = 10):
        datasets = []
        v = self.cdata.get_genepair_values(g1, g2)
        for (i,b) in enumerate(v):
            if b == self.MISSING:
                continue
            datasets.append( (self.datasets[i], self.bineffects[i][b]) )
        datasets = sorted(datasets, key=operator.itemgetter(1))
        for i in range(top):
            datasets[i] = (datasets[i][0], self.to_posterior(datasets[i][1]))

        return datasets[0:top]

    """
    Given query gene(s) return a set of connected genes/edges using
    graphle algorithm
    """
    def query(self, qgenes, edge_weight = 0, node_size = 50):
        import time
        t1 = time.time()
    

        datasize = len(self.bineffects)
        edge_weight_log = numpy.log((1.0/edge_weight)-1)

        gene_degree = {}
        for gene in qgenes:
            # Get all dataset values for gene
            values = self.cdata.get_gene_values(gene)

            t2 = time.time()
            print t2-t1

            for (i,g) in enumerate(self.genes):

                if not g:
                    continue

                gene_offset = i * datasize
                logratio = numpy.float64(self.logprior)
                for (j, bins) in enumerate(self.bineffects):
                    v = values[ gene_offset + j ]
                    if v == self.MISSING: # or use zeros bins?
                        continue
                    logratio += bins[v]

                # Filter edges by edge_weight
                if edge_weight_log < logratio:
                    continue

                ratio = numpy.exp(logratio)
                posterior = 1.0/(1.0+ratio)

                if g not in qgenes:
                    if g not in gene_degree:
                        gene_degree[g] = 0.0
                    gene_degree[g] += posterior


            t3 = time.time()
            print t3-t2




        gene_sort = sorted(gene_degree.iteritems(), key=operator.itemgetter(1), reverse=True)
        genes = set()
        genes |= qgenes
        for i in xrange(len(gene_sort[:node_size])):
            genes.add(gene_sort[i][0])

        edges = self.get_edges(genes, edge_weight)

        result = {}
        result['genes'] = []
        result['edges'] = []
        g_dict = {}
        for (i,g) in enumerate(genes):
            result['genes'].append({'query': g in qgenes, 'id' : g})
            g_dict[g] = i
        for (g1, g2, w) in edges:
            result['edges'].append({'source' : g_dict[g1], 'target' : g_dict[g2], \
                'weight' : w})

        t4 = time.time()
        print t4-t3
        return result
