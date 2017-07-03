dna = input('Enter your DNA sequence: ')
seqsize = len(dna)
if seqsize == 0:
    print('You must enter something!')
elif 0 < seqsize < 10:
    print('Your primer must have at least ten nucleotides')
elif seqsize < 25:
    print('This size is OK')
else:
    print('Your primer is too long')
