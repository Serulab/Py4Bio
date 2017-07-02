import xlrd
iedb = {}
book = xlrd.open_workbook('../../samples/sampledata.xlsx')
sh = book.sheet_by_index(0)
for row_index in range(1, sh.nrows): #skips fist line.
    iedb[int(sh.cell_value(rowx=row_index, colx=0))] = \
         sh.cell_value(rowx=row_index, colx=2)
print(iedb)
