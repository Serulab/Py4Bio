from Bio.SeqUtils import MeltingTemp as MT
import xlwt

PRIMER_FILE = '../../samples/primers.txt'
# w is the name of a newly created workbook.
w = pyExcelerator.Workbook()
# ws is the name of a new sheet in this workbook.
ws = w.add_sheet('Result')
# These two lines writes the titles of the columns.
ws.write(0, 0, 'Primer Sequence')
ws.write(0, 1, 'Tm')
i = 1
for line in open(PRIMER_FILE):
    # For each line in the input file, write the primer
    # sequence and the Tm
    j = 0
    primer = line.replace('\n','')
    ws.write(i, j, primer)
    ws.write(i, j+1, '%2.2f' %(MT.Tm_staluc(primer)))
    i += 1
# Save the spreadsheel into a file.
w.save('primerout.xls')
