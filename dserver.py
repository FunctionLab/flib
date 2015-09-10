import socket
import struct

import logging
logger = logging.getLogger(__name__)

import operator
import sys

class DataServer:

    SEARCH, QUERY = range(2)

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

    def search(self, cut, exp, didx, genes = []):
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

        scores = []
        result = s.recv(4)
        res_len = struct.unpack('<i', result)[0]

        result = s.recv(4)
        gcount = struct.unpack('<i', result)[0]

        result = s.recv(4)
        dcount = struct.unpack('<i', result)[0]

        res_len -= 8

        import time
        print 'getting results', time.time()

        # Get all bytes until finished
        str_list = []
        str_list.append( s.recv(res_len) )
        result = len(str_list[-1])

        while result < res_len:
            try:
                str_list.append( s.recv(res_len) )
                result += len(str_list[-1])
            except MemoryError:
                print len(str_list), result

        print 'starting unpack', time.time()

        try:
            bstring = ''.join(str_list)
            scores = struct.unpack('i'*dcount + 'f'*(res_len/4 - dcount), bstring)
        except AttributeError:
            print len(str_list), result
        except:
            print len(str_list), result

        print 'finished unpack', time.time()


        s.close()
        return (gcount, dcount, scores)


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
