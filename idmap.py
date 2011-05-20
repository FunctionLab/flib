import sys

class idmap:
    key_val = None
    """
    Pass the filename of the key_value pair file.
    """
    def __init__(self, filename, list=None):
        self.key_val = {}
        idfile = list
        if filename is not None:
            idfile = open(filename)

        for line in idfile:
            toks = line.strip().upper().split('\t')
            if len(toks) < 2 or toks[0] == '':
                continue

            self.key_val[toks[0]] = tuple(toks[1:])

    def keys(self):
        if self.key_val is None:
            return {}
        else:
            return self.key_val.keys()

    """
    Returns None if no file was loaded or if the key does not exist.
    """
    def get(self, id=None):
        upper_id = None
        if id is not None:
            upper_id = id.upper()
        if self.key_val is None:
            return None
        else:
            try:
                return self.key_val[upper_id]
            except KeyError:
                return None

