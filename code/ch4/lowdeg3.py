prot_seq = input('Protein sequence: ').upper()
prot_deg = {'A':4, 'C':2, 'D':2, 'E':2, 'F':2, 'G':4,
            'H':2, 'I':3, 'K':2, 'L':6, 'M':1, 'N':2,
            'P':4, 'Q':2, 'R':6, 'S':6, 'T':4, 'V':4,
            'W':1, 'Y':2}
degen_tmp = max(prot_deg.values()) * 15
for n in range(len(prot_seq) - 15):
    degen = 0
    for x in prot_seq[n:n + 15]:
        degen += prot_deg.get(x, 3.05)
    if degen <= degen_tmp:
        degen_tmp = degen
        seq = prot_seq[n:n + 15]
print(seq)
