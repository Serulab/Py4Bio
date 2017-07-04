try:
   iname = input("Enter input filename: ")
   oname = input("Enter output filename: ")
   with open(iname) as fh:
       line = fh.readline()
   if '\t' in line:
       value = line.split('\t')[0]
   with open(oname, 'w') as fw:
       fw.write(str(int(value)*.2))
except NameError:
   print("There is no TAB. Check the input file"
except FileNotFoundError:
   print("File not exist")
except PermissionError:
   print("Can't write to outfile.")
except ValueError:
   print("The value can't be converted to int")
else:
   print("Thank you!. Everything went OK.")
