import re
regex = re.compile(' |\d|\n|\t')
seq = ''
for line in open('../../samples/pMOSBlue.txt'):
    seq += regex.sub('', line)
print(seq)
