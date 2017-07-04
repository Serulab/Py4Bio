class Sequence:
    transcription_table = {'A':'U', 'T':'A', 'C':'G' , 'G':'C'}
    enz_dict = {'EcoRI':'GAATTC', 'EcoRV':'GATATC'}
    def __init__(self, seqstring):
        self.seqstring = seqstring.upper()
    def __len__(self):
        return len(self.seqstring)
    def restriction(self, enz):
        try:
            enz_target = Sequence.enz_dict[enz]
            return self.seqstring.count(enz_target)
        except KeyError:
            return 0
    def transcription(self):
        tt = ""
        for letter in self.seqstring:
            if letter in 'ATCG':
                tt += self.transcription_table[letter]
        return tt
