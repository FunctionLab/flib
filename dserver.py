import socket
import struct

import logging
logger = logging.getLogger(__name__)

import operator
import sys
import array
import numpy as np

class DataServer:

    SEARCH, QUERY, RETRIEVE = range(3)
    DSERVER_STATUS = 'dserver_status'

    def __init__(self, ip = '127.0.0.1', port = 1234):
        self.ip = ip
        self.port = port

    def open_socket(self):
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect((self.ip, self.port))
        return s

    def close_socket(self,s):
        s.shutdown(socket.SHUT_WR)
        s.close()

    def search(self, cut, exp, didx, genes = [], session = {}):

        session[self.DSERVER_STATUS] = 'Calculating results...'
        session.save()

        s = self.open_socket()

        size = 1 + 4 + 4 + 4 # opcode + dataset id + cut + exp
        size += 4*len(genes)

        size = struct.pack('<i', size)
        s.send(size)

        opcode = struct.pack('<b', self.SEARCH)
        s.send(opcode)

        did = struct.pack('<i', didx)
        s.send(did)

        params = struct.pack('<ff', cut, exp)
        s.send(params)

        gene = struct.pack('<'+'i'*len(genes), *genes)
        s.send(gene)
        s.shutdown(socket.SHUT_WR)

        result = s.recv(4)
        res_len = struct.unpack('<i', result)[0]

        result = s.recv(res_len)

        logger.debug('Total received: %s' % len(result))
        logger.debug('Total expected: %s' % (res_len))

        while len(result) < res_len:
            #logger.debug(len(result))
            result += s.recv(res_len)

        gtotal, dtotal = struct.unpack('<ii', result[0:8])

        logger.debug('Total genes %s datasets %s' % (gtotal, dtotal))

        session[self.DSERVER_STATUS] = 'Gathering results...'
        session.save()

        logger.debug('Starting unpack')

        scores = struct.unpack('<'+'f'*(dtotal + gtotal), result[8:])

        logger.debug('Finished unpack')

        s.close()

        session[self.DSERVER_STATUS] = 'Returning results...'
        session.save()

        return (gtotal, dtotal, scores)

    def retrieve(self, didx = [], genes = []):

        s = self.open_socket()

        size = 1 + 4 # opcode + dataset id + total
        size += 4*len(didx)
        size += 4*len(genes)

        size = struct.pack('<i', size)
        s.send(size)

        opcode = struct.pack('<b', self.RETRIEVE)
        s.send(opcode)

        total = struct.pack('<i', len(didx))
        s.send(total)

        did = struct.pack('<'+'i'*len(didx), *didx)
        s.send(did)

        gene = struct.pack('<'+'i'*len(genes), *genes)
        s.send(gene)
        s.shutdown(socket.SHUT_WR)

        scores = []
        result = s.recv(4)
        res_len = struct.unpack('<i', result)[0]

        # Get all bytes until finished
        result = s.recv(res_len)

        while len(result) < res_len:
            result += s.recv(res_len)

        gtotal, dtotal = struct.unpack('<ii', result[0:8])
        scores = struct.unpack('<'+'f'*gtotal*dtotal, result[8:])

        return (gtotal, scores)



if __name__ == '__main__':
    from optparse import OptionParser

    usage = "usage: %prog [options]"
    parser = OptionParser(usage, version="%prog dev-unreleased")
    parser.add_option("-I", "--IP-address",dest="ip", default='127.0.0.1', help="IP address of BNServer instance")
    parser.add_option("-p", "--port", dest="port", default=1234, help="Port number of BNServer instance", type=int)
    parser.add_option("-x", "--dataset-id", dest="did", default=0, help="Dataset ID", type=int)
    parser.add_option("-d", "--datasets", dest="dset", help="File of dataset names", metavar="FILE")
    parser.add_option("-g", "--gene-file", dest="gene_file", help="File of gene names", metavar="FILE")
    parser.add_option("-q", "--gene-query-file", dest="query_file", help="File of gene names", metavar="FILE")


    (options, args) = parser.parse_args()

    genef = open(options.gene_file)
    gidx = []
    gidx_dict = {}
    for l in genef:
        (idx, gene) = l.strip().split()
        gidx.append((int(idx)-1, gene))
        gidx_dict[gene] = int(idx) - 1
    genef.close()

    dsf = open(options.dset)
    didx = []
    for l in dsf:
        (idx, ds, pfm) = l.strip().split()
        didx.append((int(idx)-1, ds))
    dsf.close()

    qf = open(options.query_file)
    query = set()
    query_names = set()
    for l in qf:
        if l.strip() in gidx_dict:
            query.add(gidx_dict[l.strip()])
            query_names.add(l.strip())

    print query
    ds  = DataServer(options.ip, options.port)
    scores = ds.search(.5, 8, options.did, list(query))

    print len(scores)
    print scores[0:10]

#    for ((idx,name),score) in zip(gidx,genes)[0:10]:
#        print name + '\t' + ('1' if name in query_names else '-1') + '\t' + str(score)

#    for ((idx,name),score) in zip(didx,dsets)[0:10]:
#        print name + '\t' + ('1' if name in query_names else '-1') + '\t' + str(score)
