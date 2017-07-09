from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import ProtParamData
from Bio import SeqIO
with open('../../samples/pdbaa') as fh:
   for rec in SeqIO.parse(fh,'fasta'):
       myprot = ProteinAnalysis(str(rec.seq))
       print(myprot.count_amino_acids())
       print(myprot.get_amino_acids_percent())
       print(myprot.molecular_weight())
       print(myprot.aromaticity())
       print(myprot.instability_index())
       print(myprot.flexibility())
       print(myprot.isoelectric_point())
       print(myprot.secondary_structure_fraction())
       print(myprot.protein_scale(ProtParamData.kd, 9, .4))
