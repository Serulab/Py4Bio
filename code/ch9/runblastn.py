from Bio.Blast.Applications imp
BLAST_EXE = '~/opt/ncbi-blast-2
f_in = '../../samples/seq3.txt'
b_db = 'db/samples/TAIR8cds'
blastn_cline = blastn(cmd=BLAST
                      evalue=.0
rh, eh = blastn_cline()
