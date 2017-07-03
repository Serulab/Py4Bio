total_len = 0
with open('../../sampoles/B1.csv') as fh:
    next(fh)
    for n, line in enumerate(fh):
        data = line.split(',')
        total_len += int(data[1])
print(total_len/n)
