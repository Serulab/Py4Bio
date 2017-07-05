from Bio.Data.IUPACData import protein_weights as pw
protseq = input('Enter your protein sequence: ')
total_w = 0
for aa in protseq:
    total_w += pw.get(aa.upper(),0)
total_w -= 18*(len(protseq)-1)
print('The net weight is: {0}'.format(total_w))
