with open('samples/seqA.fas') as fh:
    my_file = fh.read()
name = my_file.split('\n')[0][1:]
sequence = ''.join(my_file.split('\n')[1:])
print('The name is : {0}'.format(name))
print('The sequence is: {0}'.format(sequence))
