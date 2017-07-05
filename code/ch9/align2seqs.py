from Bio.Alphabet import generic_protein
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
seq1 = 'MHQAIFIYQIGYPLKSGYIQSIRSPEYDNW'
seq2 = 'MH--IFIYQIGYALKSGYIQSIRSPEY-NW'
seq_rec_1 = SeqRecord(Seq(seq1, generic_protein), id = 'asp')
seq_rec_2 = SeqRecord(Seq(seq2, generic_protein), id = 'unk')
align = MultipleSeqAlignment([seq_rec_1, seq_rec_2])
print(align)
