import csv
data = list(csv.reader(open('../../samples/B1.csv')))
total_len = sum([int(x[1]) for x in data[1:]]))
print(total_len / (len(data)-1))
