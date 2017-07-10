import re
seq = "ATATAAGATGCGCGCGCTTATGCGCGCA"
rgx = re.compile("TAT")
i = 1
for mo in rgx.finditer(seq):
    print('Ocurrence {0}: {1}'.format(i, mo.group()))
    print('Position: From {0} to {1}'.format(mo.start(),
                                             mo.end()))
    i += 1
