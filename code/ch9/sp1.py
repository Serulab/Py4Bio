from Bio import SwissProt
with open('../../samples/spfile.txt') as fh:
    records = SwissProt.parse(fh)
    for record in records:
        print('Entry name: %s' % record.entry_name)
        print('Accession(s): %s' % ','.join(record.accessions))
        print('Keywords: %s' % ','.join(record.keywords))
        print('Sequence: %s' % record.sequence)
