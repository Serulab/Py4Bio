import csv
total_len = 0
lines = csv.reader(open('../../samples/B1.csv'))
next(lines)
for n, line in enumerate(lines):
    total_len += int(line[1])
print(total_len / n)
