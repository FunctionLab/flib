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
        
    # Not currently working when printing an object
    def __repr__(self):
        return "Counts: " + self.name
    
    def __add__(self, other):
        sum_count = Counts()
        #sum_count.name = self.name + " + " + other.name
        sum_count.num_datasets = self.num_datasets
        sum_count.tot_negats = self.tot_negats + other.tot_negats
        sum_count.tot_posits = self.tot_posits + other.tot_posits
        for i in xrange(self.num_datasets):
            sum_count.tables.append(self.tables[i] + other.tables[i])
        return(sum_count)
    
    def __mul__(self, other):
        mul_count = Counts()
        mul_count.num_datasets = self.num_datasets
        mul_count.tot_negats = self.tot_negats * other
        mul_count.tot_posits = self.tot_posits * other
        for i in xrange(self.num_datasets):
            mul_count.tables.append(self.tables[i] * other)
        return(mul_count)
    
    def __rmul__(self, other):
        return(self.__rmul__(other))

    
    def write_counts(self, filename):
        file = open(filename, 'w')
        file.write(self.name + "\t" + str(self.num_datasets) + "\n")
        file.write(str(int(self.tot_negats)) + "\t" + str(int(self.tot_posits)) + '\n')
        for i in xrange(self.num_datasets):
            file.write(self.tables[i].name + '\n')
            file.write('\t'.join([str(int(x)) for x in self.tables[i].neg_bins]) + '\n')
            file.write('\t'.join([str(int(x)) for x in self.tables[i].pos_bins]) + '\n')
        file.close()


class ConTable:
    
    def __init__(self, name="", neg_bins=None, pos_bins=None):
        self.name = name
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

