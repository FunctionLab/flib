import logging
logger = logging.getLogger(__name__)

import array
import struct

class CDatabaselet:
    def __init__(self, cdb_file, hsize, gtotal, dcount, gcount, nibble=True):
        self.header_size = hsize
        self.dataset_count = dcount
        self.gene_count = gcount
        self.cdb_file = cdb_file
        self.gene_total = gtotal

        self.nibble = nibble

        self.genes = []
        self.genes_idx = {}

    def append_gene(self,gene):
        self.genes_idx[gene] = len(self.genes)
        self.genes.append(gene)


    def __repr__(self):
        return ':'.join([self.cdb_file, str(self.header_size),
                    str(self.dataset_count), str(self.gene_count)])

    """
    Return the number of bytes for storing all datasets for a gene pair
    """
    def get_dataset_size(self):
        if self.nibble:
            return (self.dataset_count+1)/2
        else:
            return self.dataset_count

    """
    Get offset in the databaselet for gene
    """
    def get_gene_offset(self, gene):
        # header size + gene offset * # of datasets
        return self.header_size + self.genes_idx[gene] * self.get_dataset_size() * self.gene_total

    """
    Return a list of all pairwise dataset values for g1 and g2
    """
    def get_genepair_values(self, g1, g2idx):
        db_file = open(self.cdb_file, 'rb')
        seek = self.get_gene_offset(g1) + g2idx * self.get_dataset_size()
        db_file.seek(int(seek))

        byte_list = array.array('B')
        if self.nibble:
            byte_list.fromfile(db_file, self.dataset_count/2)
        else:
            byte_list.fromfile(db_file, self.dataset_count)

        if self.nibble:

            values = [None]*(len(byte_list)*2)
            for (i,b) in enumerate(byte_list):
                values[i*2] = (b & 0x0F)
                values[i*2+1] = (b >> 4)

            # check if on byte interval
            if self.nibble and int(seek) >= seek and self.dataset_count % 2 == 1:
                b = struct.unpack('B',db_file.read(1))
                values.append(b[0] & 0x0F)
        else:
            values = byte_list
        return values

    """
    Return a list of all pairwise values for g1 to all other genes across all other datasets.
    """
    def get_gene_values(self, g1):
        db_file = open(self.cdb_file, 'rb')
        seek = self.get_gene_offset(g1)
        db_file.seek(int(seek))

        byte_list = array.array('B')
        if self.nibble:
            byte_list.fromfile(db_file, (self.dataset_count + 1) * self.gene_total/2)
        else:
            byte_list.fromfile(db_file, self.dataset_count * self.gene_total)

        values = []
        if self.nibble:
            cd = 0 #cur dataset
            for b in byte_list:
                values.append( b & 0x0F )
                if cd + 1 < self.dataset_count:
                    values.append( b >> 4 )
                if cd + 2 > self.dataset_count:
                    cd = 0
                else:
                    cd += 2
        else:
            values = byte_list

        return values
