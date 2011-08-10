class Counts:

    def __init__(self):
        self.name = ""
        self.filename = ""
        self.tables = []
        self.num_datasets = 0
        self.tot_posits = 0
        self.tot_negats = 0

    def read(self, filename):
        self.filename = filename
        f = open(filename, 'r')
        line = f.readline()
        (name, num) = line.split('\t')
        self.name = name
        self.num_datasets = int(num)
        line = f.readline()
        (neg, pos) = line.split('\t')
        self.tot_negats = int(neg)
        self.tot_posits = int(pos)
        # Begin reading contingency tables
        dataset_name = f.readline().strip()
        while dataset_name:
            neg_bins = [float(x) for x in f.readline().strip().split('\t')]
            pos_bins = [float(x) for x in f.readline().strip().split('\t')]
            self.tables.append(ConTable(dataset_name, neg_bins, pos_bins))
            dataset_name = f.readline().strip()
        f.close()

    def __repr__(self):
        return "Counts: " + self.name

    def __add__(self, other):
        sum_count = Counts()
        # sum_count.name = self.name + " + " + other.name
        sum_count.num_datasets = self.num_datasets
        sum_count.tot_negats = self.tot_negats + other.tot_negats
        sum_count.tot_posits = self.tot_posits + other.tot_posits
        for i in xrange(self.num_datasets):
            sum_count.tables.append(self.tables[i] + other.tables[i])
        return(sum_count)

    def __mul__(self, other):
        mul_count = Counts()
        mul_count.name = self.name
        mul_count.num_datasets = self.num_datasets
        mul_count.tot_negats = self.tot_negats * other
        mul_count.tot_posits = self.tot_posits * other
        for i in xrange(self.num_datasets):
            mul_count.tables.append(self.tables[i] * other)
        return(mul_count)

    def __rmul__(self, other):
        return(self.__mul__(other))

    def __div__(self, other):
        return(self * (1.0 / other))

    #def __rdiv__(self, other):
        #return(self.__div__(other))

    def write_counts(self, filename, mult=1):
        file = open(filename, 'w')
        file.write(self.name + "\t" + str(self.num_datasets) + "\n")
        file.write(str(int(self.tot_negats)) + "\t" + str(int(self.tot_posits)) + '\n')
        if mult == 1:
            for i in xrange(self.num_datasets):
                file.write(self.tables[i].name + '\n')
                file.write('\t'.join([str(int(x)) for x in self.tables[i].neg_bins]) + '\n')
                file.write('\t'.join([str(int(x)) for x in self.tables[i].pos_bins]) + '\n')
        else:
            for i in xrange(self.num_datasets):
                file.write(self.tables[i].name + '\n')
                file.write('\t'.join([str(int(x * mult)) for x in self.tables[i].neg_bins]) + '\n')
                file.write('\t'.join([str(int(x * mult)) for x in self.tables[i].pos_bins]) + '\n')
        file.close()

    def to_props(self):
        for x in self.tables:
            x.to_props()

    @staticmethod
    def ave_props(counts_list, weight_list=None):
        if not len(counts_list):
            return None
        ave_counts = Counts()
        ave_counts.num_datasets = counts_list[0].num_datasets
        ave_counts.tot_negats = 1
        ave_counts.tot_posits = 1
        if weight_list:
            weight_list = [float(x) / sum(weight_list) for x in weight_list]
        for i in xrange(ave_counts.num_datasets):
            ave_counts.tables.append(ConTable.ave_props([x.tables[i] for x in counts_list], weight_list))
        return(ave_counts)

    def counts_from_props(self, mult=100000):
        self.tot_negats = mult
        self.tot_posits = mult
        for c_table in self.tables:
            c_table.neg_bins = [int(x * mult) for x in c_table.neg_props]
            c_table.pos_bins = [int(x * mult) for x in c_table.pos_props]



class ConTable:

    def __init__(self, name="", neg_bins=None, pos_bins=None):
        self.name = name
        self.neg_props = []
        self.pos_props = []
        if neg_bins:
            self.neg_bins = neg_bins[:]
        else:
            self.neg_bins = []
        if pos_bins:
            self.pos_bins = pos_bins[:]
        else:
            self.pos_bins = []

    def __add__(self, other):
        if isinstance(other, ConTable):
            sum_table = ConTable(self.name)
            for j in xrange(len(self.neg_bins)):
                sum_table.neg_bins.append(self.neg_bins[j] + other.neg_bins[j])
                sum_table.pos_bins.append(self.pos_bins[j] + other.pos_bins[j])
            return(sum_table)
        elif isinstance(other, (long, int)):
            #code for adding integers
            raise TypeError
        else:
            raise TypeError

    def __radd__(self, other):
        return(self.__add__(other))

    def to_props(self):
        tot_negats = sum(self.neg_bins)
        tot_posits = sum(self.pos_bins)
        if tot_negats == 0:
            bins = len(self.neg_bins)
            self.neg_props = [1.0 / bins] * bins
        else:
            self.neg_props = [float(x) / tot_negats for x in self.neg_bins]
        if tot_posits == 0:
            bins = len(self.neg_bins)
            self.pos_props = [1.0 / bins] * bins
        else:
            self.pos_props = [float(x) / tot_posits for x in self.pos_bins]

    @staticmethod
    def ave_props(con_table_list, weight_list=None):
        num_of_con_tables = len(con_table_list)
        num_of_bins = len(con_table_list[0].neg_props)
        c = ConTable()
        c.name = con_table_list[0].name
        if not weight_list:
            for i in range(num_of_bins):
                c.neg_props.append(sum([x.neg_props[i] for x in con_table_list]) /
                                   num_of_con_tables)
                c.pos_props.append(sum([x.pos_props[i] for x in con_table_list]) /
                                   num_of_con_tables)
        else:
            for i in range(num_of_bins):
                c.neg_props.append(sum([con_table_list[j].neg_props[i] * weight_list[j] for j in range(len(con_table_list))]))
                c.pos_props.append(sum([con_table_list[j].pos_props[i] * weight_list[j] for j in range(len(con_table_list))]))
        return(c)

    
    def __mul__(self, other):
        if isinstance(other, (long, int, float)):
            mul_table = ConTable(self.name)
            for j in xrange(len(self.neg_bins)):
                mul_table.neg_bins.append(self.neg_bins[j] * other)
                mul_table.pos_bins.append(self.pos_bins[j] * other)
            return(mul_table)
        else:
            raise TypeError

    def __rmul__(self, other):
        return(self.__mul__(other))

    def __repr__(self):
        return(self.name + '\n' +
               '\t'.join([str(x) for x in self.neg_bins]) + '\n' +
               '\t'.join([str(x) for x in self.pos_bins]) + '\n')
