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
python SplitGFF3.py --gff genome.all.gff3
Split GFF3 by sequence. The resulting file will be gff3 for each sequence, including fasta sequence at the end of gff3.
    '''
    print message

def fasta_sequence(fastafile):
    fasta_seq = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        #print record.id
        #print record.seq
        #fasta_seq[record.id] = str(record.seq)
        fasta_seq[record.id]  = record.format('fasta')
    return fasta_seq

def splitgff(gff, output):
    scaffold = defaultdict(list)
    s = re.compile(r'^##FASTA')
    s1 = re.compile(r'^#')
    flag = 0
    fastafh = open ('temp.gff3.fasta', 'w')
    with open (gff,'r') as gfffh:
        for line in gfffh:
            line = line.rstrip()
            if s.search(line):
                flag = 1
                continue
            if flag == 1:
                print >> fastafh, line
            else:
                if s1.search(line):
                    continue
                else:
                    unit = re.split(r'\t',line)
                    scaffold[unit[0]].append(line)
    fastafh.close()
    fasta_seq = fasta_sequence('temp.gff3.fasta')
    for scaf in scaffold.keys():
        outfile = output + '/' + scaf + '.gff3'
        filefh = open (outfile, 'w')
        gff_record = '\n'.join(scaffold[scaf])
        print >> filefh, '##gff-version 3'
        print >> filefh, gff_record
        print >> filefh, '##FASTA'
        #print >> filefh, '>' + scaf
        print >> filefh, fasta_seq[scaf]
        filefh.close()
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--gff')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.gff) > 0
    except:
        usage()
        sys.exit(2)

    if args.output is None:
        args.output = 'GFF3_per_sequence'
    
    os.system('mkdir ' + args.output)
    splitgff(args.gff, args.output)
    os.system('tar -zcvf ' + args.output + '.tar.gz ' + args.output)

if __name__ == '__main__':
    main()

