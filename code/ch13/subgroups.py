import re
rgx = re.compile("(?P<TBX>TATA..).*(?P<CGislands>(?:GC){3,})")
seq = "ATATAAGATGCGCGCGCTTATGCGCGCA"
result = rgx.search(seq)
print(result.group('CGislands'))
print(result.group('TBX'))
