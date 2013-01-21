import socket
import struct

import logging
logger = logging.getLogger(__name__)

import operator
import sys
import numpy
from counter import Counter

class BNServer:

    INFERENCE, DATA, GRAPH, CONTEXTS = range(4)

    def __init__(self, gidx = None, ip = '127.0.0.1', port = 1234):
        self.ip = ip
        self.port = port
        self.gidx = gidx


    def open_socket(self):
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect((self.ip, self.port))
        return s

    def close_socket(self,s):
        s.shutdown(socket.SHUT_WR)
        s.close()

    def query_network(self, genes, context = 0):
        results = {}

        s = self.open_socket()

        size = struct.pack('<i', (len(genes)+1)*4 + 1)
        s.send(size)

        opcode = struct.pack('<b', 0)
        s.send(opcode)

        contextid = struct.pack('<i', context)
        s.send(contextid)

        gene = struct.pack('<'+'i'*len(genes), *genes)
        s.send(gene)
        s.shutdown(socket.SHUT_WR)

        for g in genes:
            result = s.recv(4)
            res_len = struct.unpack('<i', result)[0]
            print res_len


            result = None
            result = s.recv(res_len)
            while len(result) < res_len:
                result += s.recv(res_len)
            res_list = struct.unpack('f'*(res_len/4), result)

            results[str(g)] = list(res_list[0:5])
            #print len(res_list), len(res_list)/float(len(self.gidx))

        s.close()
        return results

if __name__ == '__main__':
    from optparse import OptionParser

    usage = "usage: %prog [options]"
    parser = OptionParser(usage, version="%prog dev-unreleased")
    parser.add_option("-g", "--gene-file", dest="gene_file", help="File of gene names", metavar="FILE")
    parser.add_option("-i", "--IP-address",dest="ip", default='127.0.0.1', help="IP address of BNServer instance")
    parser.add_option("-p", "--port", dest="port", default=1234, help="Port number of BNServer instance", type=int)

    (options, args) = parser.parse_args()

    genef = open(options.gene_file)
    gidx = []
    for l in genef:
        (idx, gene) = l.strip().split()
        gidx.append((int(idx)-1, gene))
    genef.close()

    bns = BNServer(gidx, options.ip, options.port )
    bns.query_network([1,5,3,4], 0)
