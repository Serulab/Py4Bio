def charge_and_prop(aa_seq):
   """ Returns the net charge of a protein sequence
   and proportion of charged amino acids
   """
   protseq = aa_seq.upper()
   charge = -0.002
   cp = 0
   aa_charge = {'C':-.045, 'D':-.999, 'E':-.998, 'H':.091,
                'K':1, 'R':1, 'Y':-.001}
   for aa in protseq:
       charge += aa_charge.get(aa,0)
       if aa in aa_charge:
           cp += 1
   prop = 100.*cp/len(aa_seq)
   return (charge,prop)
