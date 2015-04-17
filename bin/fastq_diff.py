#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO

def usage():
    test="name"
    message='''
python fastq_diff.py file1.fq file2.fq

    '''
    print message

def parse_fq(fqfile):
    fqid = defaultdict(str)
    s = re.compile(r'(.*?):.*')
    for record in SeqIO.parse(fqfile,"fastq"):
        m = s.search(str(record.id))
        ids = str(record.id)
        if m:
            ids = m.groups(0)[0]
        fqid[ids] = 1
    return fqid


def readtable(infile):
    data = defaultdict(str)
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                if not data.has_key(unit[0]):
                    data[unit[0]] = unit[1]
    return data


def main():
    try:
        os.path.isfile(sys.argv[1]) and os.path.isfile(sys.argv[2])
    except:
        usage()
        sys.exit(2)

    id1 = parse_fq(sys.argv[1])
    id2 = parse_fq(sys.argv[2])
    for i in sorted(id1.keys()):
        if not id2.has_key(i):
            print i

if __name__ == '__main__':
    main()

