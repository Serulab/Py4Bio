prot_seq = input('Protein sequence: ').upper()
prot_deg = {'A':4, 'C':2, 'D':2, 'E':2, 'F':2, 'G':4,
            'H':2, 'I':3, 'K':2, 'L':6, 'M':1, 'N':2,
            'P':4, 'Q':2, 'R':6, 'S':6, 'T':4, 'V':4,
            'W':1, 'Y':2}
segs_values = []
segs_seqs = []
segment = prot_seq[:15]
a = 0
while len(segment)==15:
    degen = 0
    for x in segment:
        degen += prot_deg.get(x, 3.05)
    segs_values.append(degen)
    segs_seqs.append(segment)
    a += 1
    segment = prot_seq[a:a+15]
print(segs_seqs[segs_values.index(min(segs_values))])
