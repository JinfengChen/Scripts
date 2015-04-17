#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
import numpy as np
import re
import os
import argparse
from Bio import SeqIO
import math
from scipy import stats


def usage():
    test="name"
    message='''
python CircosConf.py --input circos.config --output pipe.conf

    '''
    print message

def fasta_id(fastafile):
    fastaid = defaultdict(str)
    for record in SeqIO.parse(fastafile,"fasta"):
        fastaid[record.id] = 1
    return fastaid


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
'''
sem: standard error of mean
two sem is the 95% confidential interval
http://www.randalolson.com/2012/08/06/statistical-analysis-made-easy-in-python/
'''
def ci95_sem(x):
    mean = np.mean(x)
    sem = stats.sem(x)
    left = mean - 2*sem
    right= mean + 2*sem
    print left, right

def ci95(x):
    mean = np.mean(x)
    sd   = np.std(x) 
    n    = len(x)
    error= 1.96*sd/math.sqrt(n)
    left = mean - error
    right= mean + error
    print left, right

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    x = [2,3,5,6,9]
    ci95(x)
    ci95_sem(x)
if __name__ == '__main__':
    main()

