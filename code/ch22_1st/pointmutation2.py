#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 19:17:25 2017

@author: sbassi
"""

from argparse import ArgumentParser

from Bio import Seq
from Bio.Alphabet import IUPAC
from Bio import Restriction
from Bio.Data import CodonTable
import jinja2

helpstr = '''Given a DNA sequence of a polypeptide, this program
generates alternative DNA sequences that code for the same 
polypeptide but that can be sorted out by DNA restriction.
Author: Sebastian Bassi (sbassi@genesdigitales.com)
License: GPL 3.0 (http://www.gnu.org/licenses/gpl-3.0.txt)'''

usage = helpstr + '\n\nusage: %(prog)s input_sequence [options]'

parser = ArgumentParser(usage=usage)
parser.add_argument('input', help='Input sequence')
parser.add_argument('-o', '--output', help='name of the output file', 
                    default='output.html')
parser.add_argument("-m", '--mutations', type=int, 
                  help='number of allowed mutations',
                  dest='n_mut', default=1)


parser.add_argument("-b", '--alignments', default=None,
                    type=int, 
                    help="alignments keep in output file")


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
                yield (letter+prox)
    
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
                try :
                    amino = t.forward_table[codon]
                except KeyError :
                    assert codon in t.stop_codons
                    continue
                try :
                    bt[amino].append(codon)
                except KeyError :
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
bakpeps = backtrans(ori_pep)

# Make a restriction analysis for the orignal sequence.
anal = Restriction.Analysis(Restriction.CommOnly, dna)
anal.print_as("map")

original_map = anal.print_that(print_=False)

# Store the enzymes that cut in the original sequence.
enzORI = anal.with_sites().keys()
enzORIset = set(enzORI)
# Get a string out of the enzyme list, only for 
# printing purposes.
oname = str(enzORI)[1:-1]
# Note: str(enzORI)[1:-1] == ", ".join(str(n) for n in enzORI)


bakpeps_out = []

bakpep_tmpdict = {}


for x in bakpeps:
    if x not in args.input:
        # Make a restriction analysis for each sequence.
        anal = Restriction.Analysis(Restriction.CommOnly, \
            Seq.Seq(x, IUPAC.unambiguous_dna))
        # Store the enzymes that cut in this sequence.
        enzTMP = anal.with_sites().keys()
        enzTMPset = set(enzTMP)
        # Get the number of mutations in backpep sequence.
        y = seqcomp(args.input, x)
        if enzTMPset!=enzORIset and enzORI!=None and y<=args.n_mut:
            # Get a string out of the enzyme list, only for 
            # printing purposes.
            all_pames.add(str(enzTMP)[1:-1])
            anal.print_as("map")
            bakpep_tmpdict['printouts'] = anal.print_that(print_=False)
            # o: Only in original sequences, p: proposed seq.
            bakpep_tmpdict['o'] = str(list(enzORIset.difference(enzTMPset)))[1:-1]
            bakpep_tmpdict['p'] = str(list(enzTMPset.difference(enzORIset)))[1:-1]
            bakpeps_out.append(bakpep_tmpdict)
            
