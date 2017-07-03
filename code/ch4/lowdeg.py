prot_seq = input('Protein sequence: ').upper()
prot_deg = {'A':4, 'C':2, 'D':2, 'E':2, 'F':2, 'G':4,
            'H':2, 'I':3, 'K':2, 'L':6, 'M':1, 'N':2,
            'P':4, 'Q':2, 'R':6, 'S':6, 'T':4, 'V':4,
            'W':1, 'Y':2}
segs_values = []
for aa in range(len(prot_seq)):
    segment = prot_seq[aa:aa + 15]
    degen = 0
    if len(segment)==15:
        for x in segment:
            degen += prot_deg.get(x, 3.05)
        segs_values.append(degen)
min_value = min(segs_values)
minpos = segs_values.index(min_value)
print(prot_seq[minpos:minpos + 15])
