import os
#import os.path
from tempfile import mkstemp

from bottle import route, post, run, static_file, request, view

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML, NCBIStandalone
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
    b_records = NCBIXML.parse(XMLin)
    for b_record in b_records:
        for alin in b_record.alignments:
            for hsp in alin.hsps:
                qs = hsp.query_start
                qe = hsp.query_end
                if qs > qe:
                    qe,qs = qs,qe
                if b_record.query not in bat1:
                    bat1[b_record.query] = [(qs,qe)]
                else:
                    bat1[b_record.query].append((qs,qe))
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
                mutable_seq[x[0]:x[1]] = mask*(x[1]-x[0])
            seq_rec = SeqRecord(mutable_seq,record.id,'','')
            outseqs.append(seq_rec)
        else:
            # Leave the sequence as found when its name is
            # not in the dictionary.
            outseqs.append(record)
    return outseqs

@route('/')
def index():
    return static_file('index.tpl', root='views/')

@post('/vector_result')
@view('result')
def result():
    seqs = request.forms.get('seqs')
    blast_db = request.forms.get('blastdb','customdb')
    if blast_db == 'customdb':
        db = os.path.join(DB_BASE_PATH, 'custom')
    elif blast_db == 'ncbivector':
        db = os.path.join(DB_BASE_PATH, 'vector')


    # Create a temporary file
    fasta_in_fh = NamedTemporaryFile(mode='w')
    # Write the user entered sequence into this temporary file
    fasta_in_fh.write(seqs)
    # Flush the data to disk without closing and deleting the file,
    # since that closing a temporary file also deletes it
    fasta_in_fh.flush()
    # Get the name of the temporary file
    file_in = fasta_in_fh.name
    # Run the BLAST query
    blastn_cline = NcbiblastnCommandline(cmd=BLAST_EXE,
                                         query=file_in,
                                         db=db,
                                         evalue=.00005,
                                         )
                                         #out="opuntia.xml")
    rh, eh = blastn_cline()

    '''
     blastn_cli = blastCL(cmd=self.blastn_path,
                             query=q_seq_fn,
                             subject=s_seq_fn,
                             task="blastn-short",
                             evalue=.00005,
                             outfmt=5, #m 7
                             out=blastout_fn)
    '''


    ##rh, eh = NCBIStandalone.blastall(BLAST_EXE, 'blastn', db,
    ##                                 file_in, expectation='1e-6')


    # Create contamination position and store it in a dictionary
    bat1 = create_rel(rh)
    # Reset the pointer position to the begining of the file
    fasta_in_fh.seek(0)
    # Get the sequences masked
    newseqs = maskseqs(fasta_in_fh, bat1)
    # Close and delete the temporary file
    fasta_in_fh.close()
    # Creates a new temporary file
    fasta_out_fh = NamedTemporaryFile()
    # Write the masked sequence into this temporary file
    SeqIO.write(newseqs,fasta_out_fh,'fasta')
    # Reset the pointer position to the begining of the file
    fasta_out_fh.seek(0)
    # Read the file
    finalout = fasta_out_fh.read()
    # Close and delete the temporary file
    fasta_out_fh.close()
    return {'finalout':finalout}

@route('/css/<filename>')
def css_static(filename):
    return static_file(filename, root='css/')

run(host='localhost', port=8080)
