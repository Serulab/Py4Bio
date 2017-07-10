import re
pattern = "[LIVM]{2}.RL[DE].{4}RLE"
with open('../../samples/Q5R5X8.fas') as fh:
    fh.readline() # Discard the first line.
    seq = ""
    for line in fh:
        seq += line.strip()
rgx = re.compile(pattern)
result = rgx.search(seq)
patternfound = result.group()
span = result.span()
leftpos = span[0]-10
if leftpos<0:
    leftpos = 0
print(seq[leftpos:span[0]].lower() + patternfound +
      seq[span[1]:span[1]+10].lower())
