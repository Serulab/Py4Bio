sequence = ''
charge = -0.002
aa_charge = {'C':-.045, 'D':-.999, 'E':-.998, 'H':.091,
            'K':1, 'R':1, 'Y':-.001}
with open('prot.fas') as fh:
   next(fh)
   for line in fh:
       sequence += line[:-1].upper()
for aa in sequence:
  charge += aa_charge.get(aa, 0)
with open('out.txt','w') as file_out:
   file_out.write(str(charge))
