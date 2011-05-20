#!/usr/bin/python

import sys
import array
import struct

class dat:
    def __init__(self, filename):
        self.gene_list = []
        self.gene_table = {}
        self.open_file(filename)
        self.gene_index = {}
        for i in range(len(self.gene_list)):
            self.gene_index[self.gene_list[i]] = i

    def open_file(self, filename):
        dab_file = open(filename, 'rb')

        #get number of genes
        a = array.array('I')
        a.fromfile(dab_file, 1)
        #print "#nodes="+str(a[0])
        size = a[0]

        #get gene names
        start = 4
        end = 4
        count = 0
        while count < a[0]:
            dab_file.seek(end)
            if(dab_file.read(2) == '\0\0'):
                dab_file.seek(start)

                gene = dab_file.read(end - start + 1)
                gene = gene.strip()
                gene = gene.replace('\x00', '')

                self.gene_list.append(gene)
                self.gene_table[gene] = count

                start = end + 2
                count += 1
                end += 1
            end += 1

        #get half matrix values
        total = (size * (size - 1)) / 2
        dab_file.seek(start)
        self.dat = array.array('f')
        self.dat.fromfile(dab_file, total)

        assert len(self.dat) == total
        #print "#edges="+str(total)
        #print self.gene_list[0:10]
        #print self.dat[0:10]

    def get_size(self):
        return len(self.gene_list)

    def get_gene(self, id):
        return self.gene_list[id]

    def get_value(self, gene1, gene2):
        #print gene1, gene2
        #try:
        #       id1 = self.gene_list.index(gene1)
        #       id2 = self.gene_list.index(gene2)
        #except ValueError:
        #       return None

        g1 = min(gene1, gene2)
        g2 = max(gene1, gene2)

        start = self.arith_sum((len(self.gene_list)) - g1, (len(self.gene_list) - 1)) #index of first id
        start += (g2 - g1) - 1 #index of second id
        try:
            v = self.dat[int(start)]
        except IndexError:
            print 'Error: ', start, gene1, gene2
            exit()

        return v
        #index=0
        #for i in range(0,len(self.gene_list)):
        #       for j in range(i+1,len(self.gene_list)):
        #               if( g1==i and g2==j ):
        #                       if( index != start ):
        #                               print index, start
        #                       assert index == start
        #                       return self.dat[index]
        #               index += 1



    def get_index(self, gene):
        try:
            return self.gene_index[gene]
        except KeyError:
            return None

    def arith_sum(self, x, y):
        return .5 * (y - x + 1) * (x + y)

    def print_table(self):
        cols = ['']
        cols.extend(self.gene_list)
        print "\t".join(cols)

        for i in range(0, self.get_size()):
            line = []
            line.append(self.gene_list[i])
            for j in range(0, i):
                v = self.get_value(i, j)
                line.append(str(v))
            line.append("0")
            for j in range(i + 1, self.get_size()):
                v = self.get_value(i, j)
                line.append(str(v))

            print "\t".join(line)


    def convert_mfinder(self):
        for i in range(0, self.get_size()):
            for j in range(i + 1, self.get_size()):
                v = self.get_value(i, j)
                if(v < .01):
                    print str(i) + " " + str(j)# + " 1"





#d = Dat(sys.argv[1])
#d.convert_mfinder()
#d.print_table()
#id1 = d.gene_list.index('MGI:3617846')
#id2 = d.gene_list.index('MGI:2685617')
#print d.get_value(id1,id2)

