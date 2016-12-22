# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 22:34:11 2016

@author: sbassi
"""

import os
from Bio.Align.Applications import ClustalwCommandline

base_dir = os.getcwd()

clustalw_exe = os.path.join(base_dir, 'clustalw2')
#cl = MultipleAlignCL('conglycinin.fasta')
data = os.path.join('samples','conglycinin.fasta')
cl = ClustalwCommandline(clustalw_exe, infile=data)
cl.outfile = 'cltest.aln'
print('Command line: {}'.format(cl))
cl()
