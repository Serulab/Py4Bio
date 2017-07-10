#!/usr/bin/env python

import argparse
import os
import sqlite3

from Bio import SeqIO, SeqRecord, Seq
from Bio.Align.Applications import ClustalwCommandline
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline as bn
from Bio import AlignIO

AT_DB_FILE = 'AT.db'
BLAST_EXE = '~/opt/ncbi-blast-2.6.0+/bin/blastn'
BLAST_DB = '~/opt/ncbi-blast-2.6.0+/db/TAIR10'
CLUSTALW_EXE = '../../clustalw2'

def allgaps(seq):
    """Return a list with tuples containing all gap positions
       and length. seq is a string."""
    gaps = []
    indash = False
    for i, c in enumerate(seq):
        if indash is False and c == '-':
            c_ini = i
            indash = True
            dashn = 0
        elif indash is True and c == '-':
            dashn += 1
        elif indash is True and c != '-':
            indash = False
            gaps.append((c_ini, dashn+1))
    return gaps

def iss(user_seq):
    """Infer Splicing Sites from a FASTA file full of EST
    sequences"""

    with open('forblast','w') as forblastfh:
        forblastfh.write(str(user_seq.seq))

    blastn_cline = bn(cmd=BLAST_EXE, query='forblast',
                      db=BLAST_DB, evalue='1e-10', outfmt=5,
                      num_descriptions='1',
                      num_alignments='1', out='outfile.xml')
    blastn_cline()
    b_record = NCBIXML.read(open('outfile.xml'))
    title = b_record.alignments[0].title
    sid = title[title.index(' ')+1 : title.index(' |')]
    # Polarity information of returned sequence.
    # 1 = normal, -1 = reverse.
    frame = b_record.alignments[0].hsps[0].frame[1]
    # Run the SQLite query
    conn = sqlite3.connect(AT_DB_FILE)
    c = conn.cursor()
    res_cur = c.execute('SELECT CDS, FULL_SEQ from seq '
                        'WHERE ID=?', (sid,))
    cds, full_seq = res_cur.fetchone()
    if cds=='':
        print('There is no matching CDS')
        exit()
    # Check sequence polarity.
    sidcds = '{0}-CDS'.format(sid)
    sidseq = '{0}-SEQ'.format(sid)
    if frame==1:
        seqCDS = SeqRecord.SeqRecord(Seq.Seq(cds),
                                     id = sidcds,
                                     name = '',
                                     description = '')
        fullseq = SeqRecord.SeqRecord(Seq.Seq(full_seq),
                                      id = sidseq,
                                      name='',
                                      description='')
    else:
        seqCDS = SeqRecord.SeqRecord(
            Seq.Seq(cds).reverse_complement(),
            id = sidcds, name='', description='')
        fullseq = SeqRecord.SeqRecord(
            Seq.Seq(full_seq).reverse_complement(),
            id = sidseq, name = '', description='')
    # A tuple with the user sequence and both AT sequences
    allseqs = (record, seqCDS, fullseq)
    with open('foralig.txt','w') as trifh:
        # Write the file with the three sequences
        SeqIO.write(allseqs, trifh, 'fasta')
    # Do the alignment:
    outfilename = '{0}.aln'.format(user_seq.id)
    cline = ClustalwCommandline(CLUSTALW_EXE,
                                infile = 'foralig.txt',
                                outfile = outfilename,
                                )
    cline()
    # Walk over all sequences and look for query sequence
    for seq in AlignIO.read(outfilename, 'clustal'):
        if user_seq.id in seq.id:
            seqstr = str(seq.seq)
            gaps = allgaps(seqstr.strip('-'))
            break
    print('Original sequence: {0}'.format(user_seq.id))
    print('\nBest match in AT CDS: {0}'.format(sid))
    acc = 0
    for i, gap in enumerate(gaps):
        print('Putative intron #{0}: Start at position {1}, '
              'length {2}'.format(i+1, gap[0]-acc, gap[1]))
        acc += gap[1]
    print('\n{0}'.format(seqstr.strip('-')))
    print('\nAlignment file: {0}\n'.format(outfilename))

description = 'Program to infer intron position based on ' \
    'Arabidopsis Thaliana genome'
parser = argparse.ArgumentParser(description=description)
ifh = 'Fasta formated file with sequence to search for introns'
parser.add_argument('input_file', help=ifh)
args = parser.parse_args()
seqhandle = open(args.input_file)
records = SeqIO.parse(seqhandle, 'fasta')
for record in records:
    iss(record)
