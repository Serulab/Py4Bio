import csv
import gzip
import os

from pymongo import MongoClient, TEXT

FILE_NAME = '../../samples/TAIR7_Transcripts_by_map_position.gz'
CONNECTION_STRING = os.getenv('MONGODB_CS', 'localhost:27017')

# Get a file handler of an uncompressed file:
with gzip.open(FILE_NAME, "rt", newline="") as f_unzip:
    rows = csv.reader(f_unzip, delimiter='\t')
    next(rows) # Skip the header
    # Dictionary for storing markers and associated information:
    at_d = {}
    # Load the dictionary using the data in the file:
    for row in rows:
        if row[0] in at_d:
            chromosome, left_val, right_val = at_d[row[0]]
            c7 = int(row[7])
            left = c7 if c7<int(left_val) else left_val
            c8 = int(row[8])
            right = c8 if c8>int(right_val) else right_val
            at_d[row[0]] = (int(chromosome), left, right)
        else:
            at_d[row[0]] = (int(row[5]), int(row[7]), int(row[8]))
# Make a list with dictionaries to be stored as documents in
# MongoDB
markers = []
for marker in at_d:
    markers.append({'marker_id': marker, 'chromosome':
                    at_d[marker][0], 'start': at_d[marker][1],
                    'end': at_d[marker][2]})
client = MongoClient(CONNECTION_STRING)
db = client.pr4
collection = db.markers_map4
collection.create_index([('marker_id', TEXT)])
collection.insert_many(markers)
