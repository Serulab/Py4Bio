import os
iname = input("Enter input filename: ")
oname = input("Enter output filename: ")
if os.path.exists(iname):
   with open(iname) as fh:
       line = fh.readline()
   if "\t" in line:
       value = line.split('\t')[0]
       if os.access(oname, os.W_OK) == 0:
           with open(oname, 'w') as fw:
               if value.isdigit():
                   fw.write(str(int(value)*.2))
               else:
                   print("Can't be converted to int")
       else:
           print("Output file is not writable")
   else:
       print("There is no TAB. Check the input file")
else:
   print("The file doesn't exist")
