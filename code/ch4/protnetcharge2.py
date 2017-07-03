prot_seq = input('Enter protein sequence: ').upper()
charge = -0.002
aa_charge = {'C':-.045,'D':-.999,'E':-.998,'H':.091,
             'K':1,'R':1,'Y':-.001}
for aa in prot_seq:
    charge += aa_charge.get('aa', 0)
print(charge)
nd{verbatim}
