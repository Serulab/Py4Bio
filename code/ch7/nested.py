iname = input("Enter input filename: ")
oname = input("Enter output filename: ")
try:
   with open(iname) as fh:
       line = fh.readline()
except FileNotFoundError:
   print("File not exist")
if '\t' in line:
   value = line.split('\t')[0]
try:
   with open(oname, 'w') as fw:
       fw.write(str(int(value)*.2))
except NameError:
   print("There is no TAB. Check the input file")
except PermissionError:
   print("Can't write to outfile.")
except ValueError:
   print("The value can't be converted to int")
else:
   print("Thank you!. Everything went OK.")
