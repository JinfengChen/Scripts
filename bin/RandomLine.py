#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse
import random

def usage():
    test="name"
    message='''
python RandomLine.py --input ../input/BGI.SNP.Jap.matrix.1 --proportion 2 --output test.file
Get random line from by proportion
    '''

    print message

def randomread(input, proportion, output):
    ofile = open (output, 'w')
    with open (input, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            r  = random.randint(1,100)
            if r <= int(proportion):
                print >> ofile, line
    ofile.close()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-p', '--proportion')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    randomread(args.input, args.proportion, args.output)

if __name__ == '__main__':
    main()

