from Bio import SwissProt
with open('../../samples/spfile.txt') as fh:
    record = next(SwissProt.parse(fh))
    for att in dir(record):
        if not att.startswith('__'):
            print(att, getattr(record, att))
