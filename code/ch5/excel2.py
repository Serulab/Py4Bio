import xlwt
list1 = [1, 2, 3, 4, 5]
list2 = [234, 267, 281, 301, 331]
wb = xlwt.Workbook()
ws = wb.add_sheet('First sheet')
ws.write(0, 0, 'Column A')
ws.write(0, 1, 'Column B')
i = 1
for x,y in zip(list1,list2): #Walk two list at the same time.
   ws.write(i, 0, x) # Row, Column, Data.
   ws.write(i, 1, y)
   i += 1
wb.save('mynewfile.xls')
