def protcharge(aa_seq):
    """Returns the net charge of a protein sequence"""
    protseq = aa_seq.upper()
    charge = -0.002
    aa_charge = {'C':-.045, 'D':-.999, 'E':-.998, 'H':.091,
                 'K':1, 'R':1, 'Y':-.001}
    for aa in protseq:
        charge += aa_charge.get(aa,0)
    return charge
