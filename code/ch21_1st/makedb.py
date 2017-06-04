import sqlite3
from Bio import SeqIO

SEQ_FILE = open('../../samples/TAIR10_seq_20101214_updated')
CDS_FILE = open('../../samples/TAIR10_cds_20101214_updated')
AT_DB_FILE = 'AT.db'

at_d = {}
# Get all sequences from TAIR sequences file.
for record in SeqIO.parse(SEQ_FILE, 'fasta'):
    sid = record.id
    seq = str(record.seq)
    at_d[sid] = [seq]
# Get all sequences from TAIR CDS file.
for record in SeqIO.parse(CDS_FILE, 'fasta'):
    sid = record.id
    seq = str(record.seq)
    at_d[sid].append(seq)
# Write to a CSV file only the entries of the dictionary that
# has data from both sources
conn = sqlite3.connect(AT_DB_FILE)
c = conn.cursor()
c.execute('create table seq(id, cds, full_seq)')
for seq_id in at_d:
    if len(at_d[seq_id])==2:
        # Write in this order: ID, CDS, FULL_SEQ.
        c.execute('INSERT INTO seq VALUES (?,?,?)',
                 ((seq_id, at_d[seq_id][1], at_d[seq_id][0])))
conn.commit()
conn.close()
