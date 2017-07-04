with open('myfile.csv') as fh:
    line = fh.readline()
value = line.split('\t')[0]
with open('other.txt',"w") as fw:
    fw.write(str(int(value)*.2))
