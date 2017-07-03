prot_seq = input('Enter your protein sequence: ')
prot_weight = {'A':89, 'V':117, 'L':131, 'I':131, 'P':115,
               'F':165, 'W':204, 'M':149, 'G':75, 'S':105,
               'C':121, 'T':119, 'Y':181, 'N':132, 'Q':146,
               'D':133, 'E':147, 'K':146, 'R':174, 'H':155}
total_weight = 0
for aa in prot_seq:
    total_weight = total_weight + prot_weight.get(aa.upper(), 0)
total_weight = total_weight - (18 * (len(prot_seq) - 1))
print('The net weight is: {0}'.format(total_weight))
