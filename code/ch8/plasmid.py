class Plasmid(Sequence):
    ab_res_dict = {'Tet':'ctagcat', 'Amp':'CACTACTG'}
    def __init__(self, seqstring):
        Sequence.__init__(self, seqstring)
    def ab_res(self, ab):
        if self.ab_res_dict[ab] in self.seqstring:
            return True
        return False
