from Bio.SeqUtils import MeltingTemp as MT

PRIMER_FILE = '../../samples/primers.txt'
for line in open(PRIMER_FILE):
    # prm stores the primer, without 5'- and -3'
    prm = line[3:len(line)-4].replace(' ','')
    # .2f is used to print up to decimals.
    print('{0},{1:.2f}'.format(prm, MT.Tm_staluc(prm)))
