import socket
import struct

import logging
logger = logging.getLogger(__name__)

import operator
import sys
import numpy
from counter import Counter
from cdatabase import CDatabase

class BNServer:

    INFERENCE, DATA, GRAPH, CONTEXTS = range(4)

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

    def inference(self, genes):
        results = {}

        s = self.open_socket()

        size = 1 + 4 # opcode + num. datasets
        for i in range(len(self.bin_effects)):
            bins = len(self.bin_effects[i])
            size += 4*(bins + 2) # data id + num. bins + bin log ratios
        size += 4*len(genes)

        size = struct.pack('<i', size)
        s.send(size)

        opcode = struct.pack('<b', 9)
        s.send(opcode)

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

    def query_network(self, genes, context = 0):
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

    cdb = CDatabase(options.cdb, didx, gidx, options.zeros_file, not options.byte)
    counter = Counter.fromcountfile(options.counts_file, cdb)

    bns = BNServer(gidx, options.ip, options.port, counter.get_bineffects())
    result = bns.inference([99])
    for g in result:
        for i in range(20):
            posterior = 1/(numpy.exp(result[g][i] + numpy.log(.99/.01)) + 1)
            print g, (i+1), posterior, result[g][i] 
