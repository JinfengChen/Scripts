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
python get_region_bam.py --input Ping.list --bam ERS467761.merge.bam > log 2> log2 &

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


#Chr1:6806761-6806763	Ping
#Chr1:36270511-36270513	Pong
def readtable(infile, bam, output):
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                mping = unit[0]
                print mping, bam
                subbam(mping, bam, output)

def subbam(mping, bam, output):
    reg = re.compile(r'(Chr\d+):(\d+)\.\.(\d+)')
    match = reg.search(mping)
    start = int(match.groups(0)[1])
    end   = int(match.groups(0)[2])
    chro    = match.groups(0)[0]
    region  = '%s:%s-%s' %(chro, start-1000, end+1000)
    #outdir  = './%s/%s' %(output, mping)
    #if not os.path.exists(outdir):
    #    os.mkdir(outdir)
    test_bam = '%s/%s.bam' %(output, mping)
    os.system('samtools view -hb %s %s > %s' %(bam, region, test_bam))
    os.system('samtools index %s' %(test_bam))

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-b', '--bam')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0 and len(args.bam) > 0
    except:
        usage()
        sys.exit(2)

    if not args.output:
        prefix = os.path.splitext(os.path.split(args.bam)[-1])[0]
        args.output = prefix
        if not os.path.exists(args.output):
            #prefix = os.path.splitext(os.path.split(args.bam)[-1])[0]
            #args.output = prefix
            os.mkdir(args.output)
    else:
        if not os.path.exists(args.output):
            os.mkdir(args.output)

    readtable(args.input, args.bam, args.output)

if __name__ == '__main__':
    main()

