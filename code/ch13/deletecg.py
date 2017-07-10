import re
regex = re.compile("(?:GC){3,}")
seq = 'ATGATCGTACTGCGCGCTTCATGTGATGCGCGCGCGCAGACTATAAG'
print('Before: {0}'.format(seq))
print('After: {0}'.format(regex.sub('', seq)))
