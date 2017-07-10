import os, io

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline as blastn
from bottle import route, post, run, static_file, request, view
from tempfile import NamedTemporaryFile

BLAST_EXE = '/home/sbassi/opt/ncbi-blast-2.6.0+/bin/blastn'
DB_BASE_PATH = '/home/sbassi/opt/ncbi-blast-2.6.0+/db/'
MASK = 'N'

def create_rel(XMLin):
    """
    Create a dictionary that relate the sequence name
    with the region to mask
    """
    bat1 = {}
    output = io.StringIO()
    output.write(XMLin)
    output.seek(0)
    b_records = NCBIXML.parse(output)
    for b_record in b_records:
        for alin in b_record.alignments:
            for hsp in alin.hsps:
                qs = hsp.query_start
                qe = hsp.query_end
                if qs > qe:
                    qe, qs = qs, qe
                record_id = b_record.query.split(' ')[0]
                if record_id not in bat1:
                    bat1[record_id] = [(qs,qe)]
                else:
                    bat1[record_id].append((qs,qe))
    return bat1

def maskseqs(ffh, bat1):
    """
    Take a FASTA file and apply the mask using the
    positions in the dictionary
    """
    outseqs = []
    for record in SeqIO.parse(ffh, 'fasta'):
        if record.id in bat1:
            # Generate a mutable sequence object to store
            # the sequence with the "mask".
            mutable_seq = record.seq.tomutable()
            coords = bat1[record.id]
            for x in coords:
                mutable_seq[x[0]:x[1]] = MASK*(x[1]-x[0])
            seq_rec = SeqRecord(mutable_seq, record.id, '', '')
            outseqs.append(seq_rec)
        else:
            # Leave the sequence as found
            outseqs.append(record)
    return outseqs

@route('/')
def index():
    return static_file('index.html', root='views/')

@post('/vector_result')
@view('vector_result')
def result():
    seqs = request.forms.get('seqs')
    db = os.path.join(DB_BASE_PATH, 'UniVec_Core')
    if request.forms.get('vector','customdb') == 'customdb':
        db = os.path.join(DB_BASE_PATH, 'custom')
    # Create a temporary file
    with NamedTemporaryFile(mode='w') as fasta_in_fh:
        # Write the user entered sequence into this temporary file
        fasta_in_fh.write(seqs)
        # Flush data to disk without closing and deleting the file,
        # since that closing a temporary file also deletes it
        fasta_in_fh.flush()
        # Get the name of the temporary file
        file_in = fasta_in_fh.name
        # Run the BLAST query
        blastn_cline = blastn(cmd=BLAST_EXE, query=file_in, db=db,
                              evalue=.0005, outfmt=5)
        rh, eh = blastn_cline()
        # Create contamination position and store it in a dictionary
        bat1 = create_rel(rh)
        # Get the sequences masked
        newseqs = maskseqs(file_in, bat1)
    with io.StringIO() as fasta_out_fh:
        SeqIO.write(newseqs, fasta_out_fh, 'fasta')
        fasta_out_fh.seek(0)
        finalout = fasta_out_fh.read()
    return {'finalout':finalout}

@route('/css/<filename>')
def css_static(filename):
    return static_file(filename, root='css/')

run(host='localhost', port=8080)
