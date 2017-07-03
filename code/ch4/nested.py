dna = input('Enter your DNA sequence: ')
seqsize = len(dna)
if seqsize < 10:
    print('Your primer must have at least ten nucleotides')
    if seqsize == 0:
        print('You must enter something!')
elif seqsize < 25:
    print('This size is OK')
else:
     print('Your primer is too long')
