import logging
logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel(logging.ERROR)

class DCheck:
    """A class for interfacing with results from DChecker from the sleipnir
    library for functional genomics.  This class exposes the properties
    positives, negatives, AUC, fpr_tpr, and rec_prec among others.
    """
    def __init__(self, dcheck=None):
        """Builds a new DCheck instance.
        
        Optional Keyword Arguments:
        dcheck -- the location of a dcheck file
        """
        self.__filename = None
        if dcheck is not None:
            self.read(dcheck=dcheck)

    @property
    def filename(self):
        "The name of the file used for this DCheck instance."
        return self.__filename

    @property
    def AUC(self):
        "The AUC from DChecker."
        return self.__AUC

    @property
    def positives(self):
        "Accessor incase we change the variable __positives."
        return self.__positives

    @property
    def negatives(self):
        "Accessor incase we change the variable __negatives."
        return self.__negatives

    @property
    def fpr_tpr(self):
        "A list of tuples of (fpr, tpr) for each cut."
        tot_pos = float(self.positives)
        tot_neg = float(self.negatives)
        return [(fp / tot_neg, tp / tot_pos) for (tp, fp, tn, fn) in self.__results]

    @property
    def rec_prec(self):
        "A list of (recall, precision) for each cut."
        tot_pos = float(self.positives)
        return [(tp / tot_pos, tp / float(tp + fp)) for (tp, fp, tn, fn) in self.__results]

    def read(self, dcheck=None):
        "Read the file located at the keyward argument dcheck and populate object."
        if dcheck is None:
            logging.error('read called with no dcheck filename.')
            return
        else:
            dfile = open(dcheck)
            for line in dfile:
                if line.startswith('#'):
                    toks = line.split('\t')
                    if toks[1] == 'P':
                        self.__positives = int(toks[2])
                    elif toks[1] == 'N':
                        self.__negatives = int(toks[2])
                    if toks[1] == 'AUC':
                        self.__results.reverse()
                        self.__AUC = float(toks[2])
                elif line.startswith('Cut'):
                    self.__results = []
                else:
                    toks = line.strip().split('\t')
                    self.__results.append(tuple([int(x) for x in toks[2:6]]))
            dfile.close()
            self.__filename = dcheck

if __name__ == '__main__':
    import sys
    from optparse import OptionParser

    usage = "usage: %prog [options]"
    parser = OptionParser(usage, version="%prog dev-unreleased")
    parser.add_option("-d", "--dcheck-file", dest="dcheck", help="file from DChecker", metavar="FILE")
    (options, args) = parser.parse_args()

    if options.dcheck is None:
        sys.stderr.write("--dcheck-file is required.\n")
        sys.exit()

    dchecker = DCheck()
    dchecker.read(options.dcheck)
    print("I just read " + dchecker.filename + ".")
    print("AUC: " + str(dchecker.AUC))
    print(dchecker.fpr_tpr)
    print(dchecker.rec_prec)
