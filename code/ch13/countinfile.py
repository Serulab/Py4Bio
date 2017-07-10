import re, sys
myregex = re.compile(sys.argv[2])
i = 0
with open(sys.argv[1]) as fh:
    for line in fh:
        i += len(myregex.findall(line))
print(i)
