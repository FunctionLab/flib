import logging
logger = logging.getLogger(__name__)

import sys
import os
import array

from cdatabaselet import CDatabaselet

class CDatabase:


    """
    Pass in location of CDatabase
    """

    def __init__(self, cdb_dir, didx, gidx, zeros_file=None, nibble=True):

        """
        Databaselets
        """
        self.cdb_genes = {}
        self.cdb_list = []

        """
        Master dataset list
        """
        self.datasets = []
        self.datasets_idx = {}
        self.datasets_zeros = []
        self.datasets_bins = []


        dblets = os.listdir(cdb_dir)
        dblets.sort()

        self.cdb = [cdb_dir + '/' + f for f in dblets]
        self.cdb_header = [None] * len(self.cdb)

        self.nibble = nibble

        """
        Setup list of datasets -- list of tuples, (idx, ds, bins)
        """
        for (idx, ds, bins) in didx:
            self.datasets_idx[str(ds)] = idx
            self.datasets.append(str(ds))
            self.datasets_bins.append(bins)

        self.datasets_zeros = [None]*len(self.datasets)

        """
        Open zeros file, if given
        """
        if zeros_file:
            zf = open(zeros_file)
            for l in zf:
                (dset, zbin) = l.strip().split('\t')
                idx = self.datasets_idx[dset]
                if idx is not None:
                    self.datasets_zeros[idx] = int(zbin)
                else:
                    logger.warning('Cannot find %s', dset)
            zf.close()

        sum_genes = 0

        """
        Load genes for each CDatabaselet
        """
        for i in range(0,len(self.cdb)):

            cdbf = open(self.cdb[i], 'rb')

            a = array.array('I')
            a.fromfile(cdbf, 4)
            self.cdb_header[i] = a[0:3]
            total_genes = a[3]
            sum_genes += total_genes

            cdbaselet = CDatabaselet(self.cdb[i], a[0], a[1], a[2], a[3], nibble)

            start, end = 16, 16
            while len(cdbaselet.genes) < total_genes:
                cdbf.seek(end)
                if( cdbf.read(1) == '\0' ):
                    cdbf.seek(start)
                    gene = cdbf.read(end-start+1).strip().replace('\x00','')
                    # store db_id and gene id in databaselet
                    self.cdb_genes[gene] = i
                    cdbaselet.append_gene(gene)
                    start = end + 1
                end += 1

            cdbf.close()
            self.cdb_list.append(cdbaselet)

        print self.datasets_zeros 
        """
        Master gene list
        """
        self.genes = [None]*sum_genes
        self.genes_idx = {}

        """
        Setup list of genes -- list of tuples, (idx, gene)
        """
        for (idx, gene) in gidx:
            self.genes_idx[str(gene)] = idx
            self.genes[idx] = str(gene)

    def is_missing(self, value):
        if self.nibble and value == 15:
            return True
        if not self.nibble and value == 255:
            return True
        else:
            return False

    def missing_val(self):
        if self.nibble:
            return 15
        else:
            return 255

    def get_dataset_zero(self, dset_idx):
        return self.datasets_zeros[dset_idx]

    def get_dataset_bins(self, dset_idx):
        return self.datasets_bins[dset_idx]

    def get_datasets(self):
        return self.datasets[:]

    def get_dataset_idx(self, dset):
        if dset in self.datasets_idx:
            return self.datasets_idx[dset]
        else:
            return None

    def gene_exists(self, gene):
        return gene in self.genes_idx

    def get_databaselet(self, gene):
        if gene in self.cdb_genes:
            return self.cdb_list[ self.cdb_genes[gene] ]
        else:
            return None

    """
    Find the databaselet for gene and return offset in the databaselet
    """

    def get_gene_offset(self, gene):
        cdata = self.get_databaselet(gene)
        return cdata.get_gene_offset(gene)

    """
    Return the byte offset in g1's databaselet at g2
    """
    def get_genes_offset(self, g1, g2):
        cdata = self.get_databaselet(g1)
        return self.get_gene_offset(g1) + self.genes_idx[g2] * cdata.get_dataset_size()

    def get_genes(self):
        return self.genes[:]

    """
    Return a list of all pairwise dataset values for g1 and g2
    """
    def get_genepair_values(self, g1, g2):
        cdata = self.get_databaselet(g1)
        if cdata is None:
            logger.warning('Could not find a CDatabaselet for gene "%s".', g1)
            return [self.missing_val()]*len(self.datasets)
        if g1 not in self.genes_idx or g2 not in self.genes_idx:
            return [self.missing_val()]*len(self.datasets)

        return cdata.get_genepair_values(g1, self.genes_idx[g2])

    """
    Return a list of all pairwise values for g1 to all other genes across all other datasets.
    """
    def get_gene_values(self, g1):
        cdata = self.get_databaselet(g1)
        if cdata is None:
            logger.warning('Could not find a CDatabaselet for gene "%s".', g1)
            return [self.missing_val()]*len(self.datasets)
        if g1 not in self.genes_idx:
            logger.warning('Could not find the gene "%s".', g1)
            return [self.missing_val()]*len(self.datasets)

        return cdata.get_gene_values(g1)

if __name__ == '__main__':
        from optparse import OptionParser

        usage = "usage: %prog [options]"
        parser = OptionParser(usage, version="%prog dev-unreleased")
        parser.add_option("-i", "--cdatabase-dir", dest="cdb", help="Directory of CDatabase", metavar="FILE")
        parser.add_option("-d", "--datasets", dest="dset", help="File of dataset names", metavar="FILE")
        parser.add_option("-g", "--gene-file", dest="gene_file", help="File of gene names", metavar="FILE")
        parser.add_option("-l", "--gene1", dest="gene1", help="Query gene")
        parser.add_option("-L", "--gene2", dest="gene2", help="Query gene")
        parser.add_option("-B", "--byte", dest="byte", help="Size of data values", action="store_true", default=False)


        (options, args) = parser.parse_args()

        if options.cdb is None:
            sys.stderr.write("--cdatabase-dir is required.\n")
            sys.exit()
        if options.gene_file is None:
            sys.stderr.write("--gene-file is required.\n")
            sys.exit()


        gene_ids = []
        for l in open(options.gene_file):
            tok = l.strip().split('\t')
            gene_ids.append((int(tok[0])-1, tok[1]))

        dset_ids = []
        for l in open(options.dset):
            tok = l.strip().split('\t')
            dset_ids.append((int(tok[0])-1, tok[1], int(tok[2])))

        #quants = None
        #if options.quant_file:
        #    quants = open(options.quant_file).readline().strip().split()
        #    print quants
        #    quants = [float(x) for x in quants]


        cdb = CDatabase(options.cdb, dset_ids, gene_ids, None, not options.byte)

        if options.gene1 and options.gene2:
            values = cdb.get_genepair_values(options.gene1, options.gene2)
            print(values)
            print(cdb.datasets)
        elif options.gene1:
            values = cdb.get_gene_values(options.gene1)
            print(values[len(cdb.datasets):2*len(cdb.datasets)])
            print(cdb.datasets)
        else:
            for g in cdb.cdb_list[0].genes:
                print(g)


