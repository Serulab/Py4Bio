dna = input('Enter your primer sequence: ')
seqsize = len(dna)
if seqsize < 10:
    print('The primer must have at least ten nucleotides')
elif seqsize < 25:
    print('This size is OK')
else:
    print('The primer is too long')
