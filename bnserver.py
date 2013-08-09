import socket
import struct

import logging
logger = logging.getLogger(__name__)

import operator
import sys
import numpy
from counter import Counter
from cdatabase import CDatabase
from xdslparser import CptNodesHolder

class BNServer:

    INFERENCE, DATA, GRAPH, CONTEXTS = range(4)
    OTF, EDGES = 9, 10


    def __init__(self, gidx = None, ip = '127.0.0.1', port = 1234, bin_effects = None):
        self.ip = ip
        self.port = port
        self.gidx = gidx
        if bin_effects:
            self.bin_effects = bin_effects


    def open_socket(self):
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect((self.ip, self.port))
        return s

    def close_socket(self,s):
        s.shutdown(socket.SHUT_WR)
        s.close()

    def get_bins_size(self):
        size = 0
        for i in range(len(self.bin_effects)):
            bins = len(self.bin_effects[i])
            size += 4*(bins + 2) # data id + num. bins + bin log ratios
        return size

    def send_bins(self, s):
        dcount = struct.pack('<i', len(self.bin_effects))
        s.send(dcount)
        for i in range(len(self.bin_effects)):
            dmessage = []
            dmessage.append(i)
            dmessage.append(len(self.bin_effects[i]))
            m = struct.pack('<'+'i'*len(dmessage), *dmessage)
            s.send(m)

            dmessage = []
            for j in range(len(self.bin_effects[i])):
                dmessage.append( self.bin_effects[i][j] )
            m = struct.pack('<'+'f'*len(dmessage), *dmessage)
            s.send(m)


    def inference_edges(self, edges):
        results = {}

        s = self.open_socket()

        size = 1 + 4 + 8 * len(edges) # opcode
        size += self.get_bins_size()

        size = struct.pack('<i', size)
        s.send(size)

        opcode = struct.pack('<b', self.EDGES)
        s.send(opcode)

        self.send_bins(s)

        size = struct.pack('<i', len(edges))
        s.send(size)

        for (g1, g2) in edges:
            e = struct.pack('<ii', g1, g2)
            s.send(e)

        s.shutdown(socket.SHUT_WR)

        result = s.recv(4)
        res_len = struct.unpack('<i', result)[0]
        result = s.recv(res_len)
        res_list = struct.unpack('f'*(res_len/4), result)

        return res_list


    def inference_otf(self, genes):
        results = {}

        s = self.open_socket()

        size = 1 + 4 # opcode + num. datasets
        #for i in range(len(self.bin_effects)):
        #    bins = len(self.bin_effects[i])
        #    size += 4*(bins + 2) # data id + num. bins + bin log ratios
        size += self.get_bins_size()
        size += 4*len(genes)

        size = struct.pack('<i', size)
        s.send(size)

        opcode = struct.pack('<b', self.OTF)
        s.send(opcode)

        self.send_bins(s)

        gene = struct.pack('<'+'i'*len(genes), *genes)
        s.send(gene)
        s.shutdown(socket.SHUT_WR)

        results = {}
        for gid in genes:
            result = s.recv(4)
            res_len = struct.unpack('<i', result)[0]

            result = None
            result = s.recv(res_len)
            while len(result) < res_len:
                result += s.recv(res_len)
            res_list = struct.unpack('f'*(res_len/4), result)
            res_list = list(res_list)
            results[gid] = res_list
            print len(results[gid])

        s.close()
        return results

    def inference(self, genes, context = 0):
        results = {}

        s = self.open_socket()

        size = struct.pack('<i', (len(genes)+1)*4 + 1)
        s.send(size)

        opcode = struct.pack('<b', self.INFERENCE)
        s.send(opcode)

        contextid = struct.pack('<i', context)
        s.send(contextid)

        gene = struct.pack('<'+'i'*len(genes), *genes)
        s.send(gene)
        s.shutdown(socket.SHUT_WR)

        for gid in genes:
            result = s.recv(4)
            res_len = struct.unpack('<i', result)[0]

            result = None
            result = s.recv(res_len)
            while len(result) < res_len:
                result += s.recv(res_len)
            res_list = struct.unpack('f'*(res_len/4), result)

            res_list = list(res_list)
            results[gid] = res_list
            #print len(res_list), len(res_list)/float(len(self.gidx))

        s.close()
        return results

    def evidence(self, gene1, gene2, prior):
        s = self.open_socket()

        size = struct.pack('<i',  9)
        s.send(size)

        opcode = struct.pack('<b', self.DATA)
        s.send(opcode)

        genes = struct.pack('<ii', gene1, gene2)
        s.send(genes)
        s.shutdown(socket.SHUT_WR)

        res_len = struct.unpack('<i', s.recv(4))[0]
        result = s.recv(res_len)
        res_list = list(struct.unpack('b'*res_len, result))

        logprior = numpy.log((1-prior)/prior)
        for (i,val) in enumerate(res_list):
            if res_list[i] != -1:
                res_list[i] = 1/(1+numpy.exp(self.bin_effects[i][val] + logprior)) - prior
            else:
                res_list[i] = None

        s.close()
        return res_list

    def data(self, gene1, gene2):
        s = self.open_socket()

        size = struct.pack('<i',  9)
        s.send(size)

        opcode = struct.pack('<b', self.DATA)
        s.send(opcode)

        genes = struct.pack('<ii', gene1, gene2)
        s.send(genes)
        s.shutdown(socket.SHUT_WR)

        res_len = struct.unpack('<i', s.recv(4))[0]
        result = s.recv(res_len)
        res_list = list(struct.unpack('b'*res_len, result))

	s.close()

        return res_list

if __name__ == '__main__':
    from optparse import OptionParser

    usage = "usage: %prog [options]"
    parser = OptionParser(usage, version="%prog dev-unreleased")
    parser.add_option("-I", "--IP-address",dest="ip", default='127.0.0.1', help="IP address of BNServer instance")
    parser.add_option("-p", "--port", dest="port", default=1234, help="Port number of BNServer instance", type=int)

    parser.add_option("-i", "--cdatabase-dir", dest="cdb", help="Directory of CDatabase", metavar="FILE")
    parser.add_option("-d", "--datasets", dest="dset", help="File of dataset names", metavar="FILE")
    parser.add_option("-g", "--gene-file", dest="gene_file", help="File of gene names", metavar="FILE")
    parser.add_option("-z", "--zeros-file", dest="zeros_file", help="File of gene names", metavar="FILE")
    parser.add_option("-f", "--counts-file", dest="counts_file", help="Counts file", metavar="FILE")
    parser.add_option("-B", "--byte", dest="byte", help="Size of data values", action="store_true", default=False)


    (options, args) = parser.parse_args()

    
    nodes = CptNodesHolder(filename=options.counts_file)

    genef = open(options.gene_file)
    gidx = []
    for l in genef:
        (idx, gene) = l.strip().split()
        gidx.append((int(idx)-1, gene))
    genef.close()

    bin_effects = []
    dsf = open(options.dset)
    didx = []
    for l in dsf:
        (idx, ds, bins) = l.strip().split()
	node = nodes.get_node(ds)
        bin_effects.append(node.get_logratios())

        didx.append((int(idx)-1, ds, int(bins)))
    dsf.close()

    #cdb = CDatabase(options.cdb, didx, gidx, options.zeros_file, not options.byte)
    #counter = Counter.fromcountfile(options.counts_file, cdb)

    bns = BNServer(gidx, options.ip, options.port, bin_effects)
    #result = bns.inference_otf([1])
    #for g in result:
    #    for i in range(1,2):
            #posterior = 1/(numpy.exp(result[g][i] + numpy.log(.99/.01)) + 1)
    #        print g, (i+1), result[g][i]

    print bns.data(2,1)
    print bns.data(1,2)
    #lr = filter( lambda x: x is not None, bns.evidence(96,1,.01) )
    #print numpy.sum(lr)
