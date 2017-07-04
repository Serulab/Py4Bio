class Sequence:
    transcription_table = {'A':'U', 'T':'A', 'C':'G', 'G':'C'}
    comp_table = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    def __init__(self, seqstring):
        self.seqstring = seqstring.upper()
    def restriction(self, enz):
        enz_dict = {'EcoRI':'ACTGG', 'EcoRV':'AGTGC'}
        try:
            target = enz_dict[enz]
        except KeyError:
            raise ValueError('No such enzime in out enzime DB')
        return self.seqstring.count(target)
    def __getitem__(self,index):
        return self.seqstring[index]
    def __getslice__(self, low, high):
        return self.seqstring[low:high]
    def __len__(self):
        return len(self.seqstring)
    def __str__(self):
        if len(self.seqstring) >= 28:
            return '{0}...{1}'.format(self.seqstring[:25],
                                      self.seqstring[-3:])
        else:
            return self.seqstring
    def transcription(self):
        tt = ''
        for x in self.seqstring:
            if x in 'ATCG':
                tt += self.transcription_table[x]
        return tt
    def complement(self):
        tt = ''
        for x in self.seqstring:
            if x in 'ATCG':
                tt += self.comp_table[x]
        return tt
