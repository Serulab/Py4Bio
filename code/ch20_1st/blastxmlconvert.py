#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
From BLAST XML to HTML. By Sebastian Bassi.
Tested with BLASTN xml files from 2.2.16 to 2.2.18.
BLASTN xml files < 2.2.16 are not properly formatted.
Converts a single BLAST XML to one HTML file.

Created on Tue Jan 17 20:30:10 2017

@author: sbassi
"""

from argparse import ArgumentParser
import xml.etree.cElementTree as cET

from jinja2 import Template

helpstr = '''XML2HTML converts a BLAST XML file into one, 
or multiple HTML files. 

Author: Sebastian Bassi (sbassi@genesdigitales.com)
Thanks to Yoan Jacquemin for help in testing.
License: GPL 3.0 (http://www.gnu.org/licenses/gpl-3.0.txt)'''

usage = helpstr + '\n\nusage: %(prog)s input_file [options]'

parser = ArgumentParser(usage=usage)
parser.add_argument('input', help='Input file name (in XML format)',
                    nargs='+')
parser.add_argument('-o', '--output', help='name of the output file', 
                    default='output.html')
parser.add_argument("-d", '--descriptions', type=int, 
                  help="descriptions keep in output file",
                  dest='desc')
parser.add_argument("-b", '--alignments', default=None,
                    type=int, 
                    help="alignments keep in output file")

def get_data(f_in):
    """
    """
    data = {}
    tree = cET.parse(f_in)
    root = tree.getroot()
    # Get head data
    data['version_date'] = root.find('BlastOutput_version').text
    data['application'] = root.find('BlastOutput_program').text
    data['reference'] = root.find('BlastOutput_reference').text[12:]
 
    
    
    return data

def htmlhead(f_in,outf):
    tree = cET.parse(f_in)
    root = tree.getroot()
    version_date = root.find('BlastOutput_version').text
    application = root.find('BlastOutput_program').text
    reference = root.find('BlastOutput_reference').text[12:]
    fo = open(outf,'w')
    fo.write('''<HTML>
<TITLE>BLAST Search Results</TITLE>
<BODY BGCOLOR="#FFFFFF" LINK="#0000FF" \
VLINK="#660099" ALINK="#660099">
<!-- Generated from %s by XML2HTML (Sebastian Bassi) -->
<PRE>''' %(f_in))
    fo.write('<b>%s %s</b>' %(application,version_date))
    fo.write('\n<b><a href="http://www.ncbi.nlm.nih.gov/entrez/\
query.fcgi?db=PubMed&cmd=Retrieve&list_uids=9254694&dopt=\
Citation">Reference</a>:</b>'+reference.replace('~',' ')
+'\n')
    return (fo,root)

def htmlfoot(fo,f_in,root):
    b_version = root.findtext('BlastOutput_db')
    num_letter_db = root.findtext(
    'BlastOutput_iterations/Iteration/Iteration_stat/Statistics/\
Statistics_db-num')
    num_seqs_db = root.findtext(
    'BlastOutput_iterations/Iteration/Iteration_stat/Statistics/\
Statistics_db-len')
    lambd = root.findtext(
    'BlastOutput_iterations/Iteration/Iteration_stat/Statistics/\
Statistics_lambda')
    kappa = root.findtext(
    'BlastOutput_iterations/Iteration/Iteration_stat/Statistics/\
Statistics_kappa')
    entrop = root.findtext(
    'BlastOutput_iterations/Iteration/Iteration_stat/Statistics/\
Statistics_entropy')
    b_prg = root.findtext('BlastOutput_program')
    p_sc_match = root.findtext('BlastOutput_param/Parameters/\
Parameters_sc-match')
    p_sc_mismatch = root.findtext('BlastOutput_param/Parameters/\
Parameters_sc-mismatch')
    p_gap_open = root.findtext('BlastOutput_param/Parameters/\
Parameters_gap-open')
    p_gap_extend = root.findtext('BlastOutput_param/Parameters/\
Parameters_gap-extend')
    fo.write('''<PRE>
  Database: %s
  Number of letters in database: %s
  Number of sequences in database:  %s

Lambda     K      H
    %.2f    %.3f     %.2f

Matrix: %s matrix:%s %s
Gap Penalties: Existence: %s, Extension: %s
Number of Sequences: %s
Length of database: %s
</PRE>
</BODY>
</HTML>''' %(b_version,num_letter_db,num_seqs_db,
             float(lambd),float(kappa),
             float(entrop),b_prg,p_sc_match,
             p_sc_mismatch,p_gap_open,p_gap_extend,
             num_seqs_db,num_letter_db))
    fo.close()
    return None

def prettyalign(fo,q,qs,qe,m,s,ss,se):
    """ Format the alignment in slices of 60 characters
    """
    #fo=file handler
    #q query sequence
    #qs query_start (or query_from)
    #qe query_end (or query_to)
    #m match sequence
    #s, ss and se are the equivalent for subject/hit
    pos = 0
    qr=range(qs,qe-1,-1) if qs>qe else range(int(qs),int(qe)+61)
    qini = qs
    qend = qe
    sr = range(ss,se-1,-1) if ss>se else range(ss,ss+len(s))
    mx = max(len(str(qr[-1])),len(str(sr[-1])))
    q_desp = 0
    s_desp = 0
    if max(len(q),len(s))>=60:
        finant_u = (pos+1 if ss>se else pos-1)
        finant_d = (pos+1 if ss>se else pos-1)
        while pos<max(len(q)-(len(q)%60),len(s)-(len(s)%60)):
            q_desp += (q[pos:pos+60].count('-')
                       if '-' in q[pos:pos+60] else 0)
            s_desp += (s[pos:pos+60].count('-')
                       if '-' in s[pos:pos+60] else 0)
            fo.write('Query: %-*s %s %s\n'%(mx,
                     qr[finant_u-1 if ss>se else finant_u+1],
                     q[pos:pos+60],qr[pos+59-q_desp]))
            fo.write('       '+' '*mx+' '+m[pos:pos+60]+'\n')
            fo.write('Sbjct: %-*s %s %s\n\n'%
                     (mx,sr[finant_d-1 if ss>se else finant_d+1],
                     s[pos:pos+60],sr[pos+59-s_desp]))
            finant_u = pos+59-q_desp
            finant_d = pos+59-s_desp
            pos += 60
        if len(q)%60!=0:
            q_desp += (q[pos:pos+60].count('-')
                       if '-' in q[pos:pos+60] else 0)
            s_desp += (s[pos:pos+60].count('-')
                       if '-' in s[pos:pos+60] else 0)
            fo.write('Query: %-*s %s %s\n'%(mx,qr[pos-q_desp],
                     q[pos:pos+60],qend))
            fo.write('       '+' '*mx+' '+m[pos:pos+60]+'\n')
            fo.write('Sbjct: %-*s %s %s\n\n'%(mx,sr[pos-s_desp],
                     s[pos:pos+60],sr[-1]))
    else:
        fo.write('Query: %-*s %s %s\n'%(mx,qini,
                 q[pos:pos+60],qend))
        fo.write('       '+' '*mx+' '+m[pos:pos+60]+'\n')
        fo.write('Sbjct: %-*s %s %s\n\n'%(mx,sr[pos],
                 s[pos:pos+60],sr[-1]))
    return None

def blastconv(f_in,fo,de,al):
    i_hits = {}
    hits = {}
    hsps = {}
    for ev,x in cET.iterparse(f_in):
        if 'BlastOutput_query-def' in x.tag:
            b_query_def = x.text
        elif 'BlastOutput_query-len' in x.tag:
            b_query_len = x.text
        elif 'BlastOutput_db' in x.tag:
            b_db = x.text
        elif 'Statistics_db-num' in x.tag:
            s_db_num=x.text
        elif 'Statistics_db-len' in x.tag:
            s_db_len=x.text
        elif 'Parameters_expect' in x.tag:
            p_expect = x.text
        elif 'Parameters_filter' in x.tag:
            p_filter = x.text
        elif 'Iteration_query-def' in x.tag:
            i_query_def = x.text
        elif 'Iteration_iter-num' in x.tag:
            i_iter_num = x.text
        elif 'Iteration_query-ID' in x.tag:
            i_query_id = x.text
        elif 'Iteration_query-len' in x.tag:
            i_query_len = x.text
        elif 'Iteration'==x.tag:
            i_hits[int(i_iter_num)] = (i_query_id, i_query_def,
                                       i_query_len, hits)
            hits = {}
        elif 'Hit_num' in x.tag:
            h_num = x.text
        elif 'Hit_id' in x.tag:
            h_id = x.text
        elif 'Hit_def' in x.tag:
            h_def = x.text
        elif 'Hit_accession' in x.tag:
            h_accession = x.text
        elif 'Hit_len' in x.tag:
            h_len = x.text
        elif 'Hit'==x.tag:
            hits[int(h_num)] = (h_id,h_def,h_accession,h_len,
                                hsps)
            hsps = {}
        elif 'Hsp_num' in x.tag:
            hsp_num = x.text
        elif 'Hsp_bit-score' in x.tag:
            hsp_bit_score = x.text
        elif 'Hsp_score' in x.tag:
            hsp_score = x.text
        elif 'Hsp_evalue' in x.tag:
            hsp_evalue = x.text
        elif 'Hsp_query-from' in x.tag:
            hsp_query_from = x.text
        elif 'Hsp_query-to' in x.tag:
            hsp_query_to = x.text
        elif 'Hsp_hit-from' in x.tag:
            hsp_hit_from = x.text
        elif 'Hsp_hit-to' in x.tag:
            hsp_hit_to = x.text
        elif 'Hsp_query-frame' in x.tag:
            hsp_query_frame = x.text
        elif 'Hsp_hit-frame' in x.tag:
            hsp_hit_frame = x.text
        elif 'Hsp_identity' in x.tag:
            hsp_identity = x.text
        elif 'Hsp_positive' in x.tag:
            hsp_positive = x.text
        elif 'Hsp_align-len' in x.tag:
            hsp_align_len = x.text
        elif 'Hsp_qseq' in x.tag:
            hsp_qseq = x.text
        elif 'Hsp_hseq' in x.tag:
            hsp_hseq = x.text
        elif 'Hsp_midline' in x.tag:
            hsp_mid = x.text
        elif 'Hsp'==x.tag:
            try:
                hn = (hsp_bit_score,hsp_score,hsp_evalue,
                      hsp_query_from,hsp_query_to,
                      hsp_hit_from,hsp_hit_to,
                      hsp_query_frame,hsp_hit_frame,
                      hsp_identity,hsp_positive,
                      hsp_align_len,hsp_qseq,hsp_hseq,hsp_mid)
            except UnboundLocalError:
                hn = (hsp_bit_score,hsp_score,hsp_evalue,
                      hsp_query_from,hsp_query_to,
                      hsp_hit_from,hsp_hit_to,
                      hsp_query_frame,hsp_query_frame,
                      hsp_identity,hsp_positive,
                      hsp_align_len,hsp_qseq,hsp_hseq,hsp_mid)
            hsps[int(hsp_num)]= hn
        elif 'Statistics_hsp-len' in x.tag:
            s_hsp_len = x.text
    ihits = list(i_hits.keys())
    ihits.sort()
    # Iterations with no BLAST result are missing, so the script
    # can't iterate over a range from one to the end.
    for x in ihits:
        fo.write('<b>Query=</b> %s\n' %(i_hits[x][1]))
        fo.write('         (%s letters)\n' %(i_hits[x][2]))
        fo.write('<b>Database:</b> %s\n' %(b_db))
        fo.write('         %s sequences; %s total letters\n'
                 %(s_db_num,s_db_len))
        fo.write('''Searching...................................\
...done
<PRE>


                                                             \
Score    E
Sequences producing significant alignments:                  \
(bits) Value

''')
        for y in range(1,len(i_hits[x][3])+1)[:de]:
            k = i_hits[x][3][y][0]
            desc = i_hits[x][3][y][1]
            bs = i_hits[x][3][y][4][1][0]
            sc = i_hits[x][3][y][4][1][2]
            if 'gi|' in k:
                m = k.index('gi|')+3
                gi = k[m:k[m:].index('|')+3]
                fo.write('<a href="http://www.ncbi.nlm.nih.gov/\
entrez/query.fcgi?cmd=Retrieve&db=Nucleotide&list_uids=%s&dopt=\
GenBank" >%s</a> %s <a href = #%s> %s</a>   %s\n' %(gi,
                k.replace('gi|'+gi+'|',''),desc[:36]+'...',gi,
                          bs,sc))
            else:
                fo.write('><a name=%s></a>%s\n'%(k,desc))
        fo.write('\n</PRE>\n')
        for y in range(1,len(i_hits[x][3])+1)[:al]:
            fo.write('<PRE>\n')
            k = i_hits[x][3][y][0]
            desc = i_hits[x][3][y][1]
            if 'gi|' in k:
                m = k.index('gi|')+3
                gi = k[m:k[m:].index('|')+3]
                fo.write('><a href="http://www.ncbi.nlm.nih.gov/\
entrez/query.fcgi?cmd=Retrieve&db=Nucleotide&list_uids=%s&dopt=\
GenBank" >%s</a> %s \n' %(gi,k.replace('gi|'+gi+'|',''),desc))
            else:
                fo.write('>%s %s\n' %(k,desc))
            fo.write(' Length = '+i_hits[x][3][y][3]+'\n')
            # Walk over all the hsps
            for z in range(1,len(i_hits[x][3][y][4])+1):
                bs = i_hits[x][3][y][4][z][0]
                hsc = i_hits[x][3][y][4][z][1]
                sc = i_hits[x][3][y][4][z][2]
                h_id = i_hits[x][3][y][4][z][10]
                h_pos = i_hits[x][3][y][4][z][11]
                q_frame = i_hits[x][3][y][4][z][7]
                h_frame = i_hits[x][3][y][4][z][8]
                q_from = int(i_hits[x][3][y][4][z][3])
                q_to = int(i_hits[x][3][y][4][z][4])
                h_from = int(i_hits[x][3][y][4][z][5])
                h_to = int(i_hits[x][3][y][4][z][6])
                qseq = i_hits[x][3][y][4][z][12]
                hseq = i_hits[x][3][y][4][z][13]
                mid = i_hits[x][3][y][4][z][14]
                qf = 'Plus' if int(q_frame)>0 else 'Minus'
                hf = 'Plus' if int(h_frame)>0 else 'Minus'
                fo.write('Score = %s bits (%s), Expect = %s\n'
                         %(bs,hsc,sc))
                fo.write('Identities = %s/%s (%.0f%%)\n'
                         %(h_id,h_pos,
                           float(int(h_id))/int(h_pos)*100))
                fo.write('Strand = %s/%s\n\n\n' %(qf,hf))
                prettyalign(fo,qseq,q_from,q_to,mid,hseq,
                            h_from,h_to)
                fo.write('</PRE>\n')
    return fo

def doconvert(f_in,outfile,desc,align):
    fo,root = htmlhead(f_in,outfile,)
    fo = blastconv(f_in,fo,desc,align)
    htmlfoot(fo,f_in,root)
    return None

args = parser.parse_args()

if len(args.input) == 1:
    f = args.input[0]
    if args.output is None:
        args.output = f[:-3]+'html'
    doconvert(f, args.output, args.desc, args.alignments)
elif len(args.input)>1:
    for f in args:
        outfile = f[:-3]+'html'
        doconvert(f, outfile, args.desc, args.alignments)