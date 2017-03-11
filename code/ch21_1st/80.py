#!/usr/bin/env python

import sys
import os
from Bio import SeqIO, SeqRecord, Seq, Clustalw
from Bio.Blast import NCBIStandalone
from Bio.Blast import NCBIXML
from Bio.Clustalw import MultipleAlignCL

AT_DB_FILE = 'AT.db'
blast_exe ='/home/sb/blast-2.2.20/bin/blastall'
blast_db = '/home/sb/blast-2.2.20/bin/TAIR8cds'

def allgaps(seq):
    """Return a list with tuples containing all gap positions
       and length. seq is a string."""
    i = 0
    gaps = []
    indash = False
    for c in seq:
        if indash is False and c=='-':
            c_ini = i
            indash = True
            dashn = 0
        elif indash is True and c=='-':
            dashn += 1
        elif indash is True and c!='-':
            indash = False
            gaps.append((c_ini,dashn+1))
        i += 1
    return gaps

def iss(record):
    """Infer Splicing Sites from a FASTA file full of EST
    sequences"""

    usersid = record.id
    userseq = record.seq
    result, err = NCBIStandalone.blastall(blast_exe, "blastn",
                  blast_db, f_name, expectation='1e-10',
                  descriptions='1', alignments='1')

    of = open('outfile.xml','w')
    of.write(result.read())
    result.close()
    of.close()
    b_record = NCBIXML.parse(open('outfile.xml')).next()
    title = b_record.alignments[0].title
    sid = title[title.index(' ')+1:title.index(' |')]

    # Polarity information of returned sequence.
    # 1 = normal, -1 = reverse.
    frame = b_record.alignments[0].hsps[0].frame[1]

    # Run the SQLite query
    ###NO!!
    conn = sqlite3.connect(AT_DB_FILE)
    c = conn.cursor()
    print(c.execute('SELECT * from seqs WHERE ID=?', sid))
    xxx

    result = x.readline().split('|')
    cds = result[1]
    seq = result[2][:-1]

    if cds=='':
        print 'There is no matching CDS'
        exit()

    # Check sequence polarity.
    if frame==1:
        seqCDS = SeqRecord.SeqRecord(Seq.Seq(cds),id=sid+'-CDS'
                                 ,name="",description="")
        fullseq = SeqRecord.SeqRecord(Seq.Seq(seq),id=sid+'-SEQ'
                                 ,name="",description="")
    else:
        seqCDS = SeqRecord.SeqRecord(
            Seq.Seq(cds).reverse_complement(),id=sid+'-CDS',
            name="",description="")
        fullseq = SeqRecord.SeqRecord(
            Seq.Seq(seq).reverse_complement(),id=sid+'-SEQ',
            name="",description="")

    # Create a tuple with the user sequence and both AT sequences.
    allseqs = (record,seqCDS,fullseq)

    trifh = open('foralig.txt','w')
    # Write the file with the three sequences.
    SeqIO.write(allseqs,trifh,"fasta")
    trifh.close()

    # Do the alignment:
    cline = MultipleAlignCL('foralig.txt')
    cline.set_output(usersid+".aln")
    alignment = Clustalw.do_alignment(cline)

    # Walk over all aligned sequences and look for query sequence
    for seq in alignment.get_all_seqs():
        if usersid in seq.id:
            seqstr = seq.seq.tostring()
            gaps = allgaps(seqstr.strip('-'))
            break

    print "Original sequence:",usersid
    print "\nBest match in AT CDS:",sid

    i = 0
    acc = 0
    for gap in gaps:
        i += 1
        print "Intron #%s: Start at position %s, length %s"\
              %(i,gap[0]-acc,gap[1])
        acc += gap[1]

    print '\n'+seqstr.strip('-')
    print '\nAlignment file: '+usersid+'.aln\n'
    return None

try:
    f_name = sys.argv[1]
except:
    print "Run this program from command line as:"
    print "iss.py file_in"
    exit()

## DEBUG: f_name='/mnt/hda2/bio/t3.txt'
seqhandle = open(f_name)
records = SeqIO.parse(seqhandle, "fasta")

for record in records:
    iss(record)
This code is part of the book "Python for Bioinformatics", by Sebastian Bassi (sbassi@genesdigitales.com). Return to home page.
