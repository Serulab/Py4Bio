sequence = ''
with open('samples/seqA.fas') as fh:
    name = fh.readline()[1:-1]
    for line in fh:
        sequence += line.replace('\n','')
print('The name is : {0}'.format(name))
print('The sequence is: {0}'.format(sequence))
