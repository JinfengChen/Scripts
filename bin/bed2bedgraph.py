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
python bed2bedgraph.py --input test.bed

bedtools bamtobed -i ~/BigData/04.Epigenome/input/mPing/Jinghua20140319/Sample_MZ19.bam -tag > H3K9me.tag.bed &
bed:
Chr1    1000    1100    HAL:1303:C30N5ACXX:6:1216:2101:11496    255     +
Chr1    1001    1101    HAL:1303:C30N5ACXX:6:1115:11392:36761   255     +
Chr1    1001    1101    HAL:1303:C30N5ACXX:6:1211:12115:76601   255     -
bedgraph:
Chr1	1000	1100	1
Chr1	1001	1101	2'''    
    print message

'Chr1    1000    1100    HAL:1303:C30N5ACXX:6:1216:2101:11496    255     +'
def bed2graph(infile):
    data = defaultdict(int)
    last = ''
    print 'track type=bedGraph name="H3K9me2"';
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            if len(line) > 2: 
                unit = re.split(r'\t',line)
                index = unit[0] + '_' + unit[1]
                if index == last:
                    data[index] += 1
                else:
                    data[index] = 1
                    print '%s\t%s\t%s\t%s' %(unit[0], unit[1], unit[2], data[last])
                    del data[last]
                    last = index

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

    bed2graph(args.input)

if __name__ == '__main__':
    main()

