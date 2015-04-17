#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
python SplitFasta.py --fasta /rhome/cjinfeng/BigData/00.RD/seqlib/MSU_r7.fa --outdir MSU7_unmask
Split genome fasta file and name every sequence as the name of sequence. Output will be in a new dir: outdir
    '''
    print message

def split(fasta, outdir):
    for record in SeqIO.parse(fasta,"fasta"):
        ofile = open (outdir + '/' + record.id + '.fa', 'w')
        SeqIO.write(record,ofile,"fasta")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta')
    parser.add_argument('-o', '--outdir')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.fasta) > 0
    except:
        usage()
        sys.exit(2)

    if args.outdir is None:
        args.outdir = 'Split'
   
    os.system('mkdir %s' % (args.outdir))
    split(args.fasta, args.outdir)

if __name__ == '__main__':
    main()

