import re, sys
myregex = re.compile(sys.argv[2])
counter = 0
with open(sys.argv[1]) as fh:
    for line in fh:
        if myregex.search(line):
            counter += 1
print(counter)
