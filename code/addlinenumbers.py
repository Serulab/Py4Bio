#!/usr/bin/env python

from __future__ import print_function
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('filein')
args = parser.parse_args()

with open(args.filein) as filein:
    file_len = (len(filein.readlines()))
padding = 2 if 99 > file_len > 9 else 1
ln = 1
if padding <=2:
    with open(args.filein) as filein:
        for line in filein:
            if ln < 10:
                print(' '*(padding+1) + str(ln) + ' ' + line, end='')
            else:
                print(' '*padding + str(ln) + ' ' + line, end='')
            ln += 1
else:
    with open(args.filein) as filein:
        for line in filein:
            print(' '*4 + str(ln) + ' ' + line, end='')
            ln += 1
