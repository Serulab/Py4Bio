from argparse import ArgumentParser

from Bio import Seq
from Bio.Alphabet import IUPAC
from Bio import Restriction
from Bio.Data import CodonTable
from jinja2 import Template

helpstr = '''Given a DNA sequence of a polypeptide, this program
generates alternative DNA sequences that code for the same
polypeptide but that can be sorted out by DNA restriction.
Requires Biopython >= 1.69
Author: Sebastian Bassi (sbassi@genesdigitales.com)'''
usage = helpstr + '\n\nusage: %(prog)s input_sequence [options]'
parser = ArgumentParser(usage=usage)
parser.add_argument('input', help='Input sequence')
parser.add_argument('-o', '--output', help=
                    'name of the output file',
                    default='output.txt')
parser.add_argument('-m', '--mutations', type=int,
                  help='number of allowed mutations',
                  dest='n_mut', default=1)
parser.add_argument('-t', '--table', type=int,
                  help='translation table',
                  dest='table_id', default=1)

def backtrans(ori_pep, table_id=1):
    """
    Function to make backtranslation (from peptide to DNA)
    This function needs the peptide sequence and the code of
    translation table. Code number is the same as posted in:
    http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
    """
    def recurs(order, pos):
        for letter in bt[order[pos]]:
            if pos == len(order) - 1:
                yield letter
                continue
            for prox in recurs(order, pos+1):
                yield (letter + prox)

    def combine(order):
        ordened = set()
        for frase in recurs(order, 0):
            ordened.add(frase)
        return ordened

    t = CodonTable.generic_by_id[table_id]
    bt = dict()
    for a1 in "ATCG" :
        for a2 in "ATCG" :
            for a3 in "ATCG" :
                codon = a1 + a2 + a3
                try:
                    amino = t.forward_table[codon]
                except KeyError:
                    assert codon in t.stop_codons
                    continue
                try:
                    bt[amino].append(codon)
                except KeyError:
                    bt[amino] = [codon]
    return list(combine(ori_pep))

def seqcomp(s1, s2):
    """
    Compares 2 sequences and returns a value with
    how many differents elements they have.
    """
    p = len(s1)
    for x,y in zip(s1, s2): # Walk through 2 sequences.
        if x==y:
            p -= 1
    return p

args = parser.parse_args()
dna = Seq.Seq(args.input, IUPAC.unambiguous_dna)
# Translate DNA sequence.
ori_pep = dna.translate()
# Get all backtranslations.
bakpeps = backtrans(ori_pep, args.table_id)
# Make a restriction analysis for the orignal sequence.
analysis = Restriction.Analysis(Restriction.CommOnly, dna)
analysis.print_as('map')
ori_map = analysis.format_output()
# Store the enzymes that cut in the original sequence.
enz = list(analysis.with_sites().keys())
# Get a string out of the enzyme list, for printing.
oname = str(enz)[1:-1]
enz = set(enz)
bakpeps_out = []
for bakpep in bakpeps:
    tmp_d = {}
    if bakpep not in args.input:
        # Make a restriction analysis for each sequence.
        analysis = Restriction.Analysis(Restriction.CommOnly,
                Seq.Seq(bakpep, IUPAC.unambiguous_dna))
        # Store the enzymes that cut in this sequence.
        enz_tmp = list(analysis.with_sites().keys())
        enz_tmp = set(enz_tmp)
        # Get the number of mutations in backpep sequence.
        y = seqcomp(args.input, bakpep)
        if enz_tmp != enz and enz and y <= args.n_mut:
            analysis.print_as('map')
            tmp_d['pames'] = str(enz_tmp)[1:-1]
            tmp_d['graph'] = analysis.format_output()
            tmp_d['ori_seq'] = str(list(enz.difference(
                                        enz_tmp)))[1:-1]
            tmp_d['proposed_seq'] = str(list(
                        enz_tmp.difference(enz)))[1:-1]
            bakpeps_out.append(tmp_d)
with open('mutation.tpl') as fh:
    template = Template(fh.read())
    render = template.render(dna_input=args.input,
                             ori_pep=ori_pep,
                             ori_map=ori_map,
                             oname=oname,
                             bakpeps_out=bakpeps_out)
with open(args.output, 'w') as file_out:
    file_out.write(render)
