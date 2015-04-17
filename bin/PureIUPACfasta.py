#!/opt/Python/2.7.3/bin/python
import sys
from collections import defaultdict
from numpy import *
import re
import os
import argparse

def usage():
    test="name"
    message='''
python PureIUPACfasta.py --input test.fa
Convert IUPAC code in fasta into ATCG
test                 prefix
test.fa              fasta
test.temp.fa         temp file that convert iupac into atcg in lower case
test.noIUPAC.fa      final file convert all base into ATCG
    '''
    print message

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.input) > 0
    except:
        usage()
        sys.exit(2)

    s = re.compile(r'(.*)\.(f\w+a)')
    m = s.search(args.input)
    if m:
        prefix = m.groups(0)[0]
        sufix  = m.groups(0)[1]
        temp   = prefix + '.temp.' + sufix
        final  = prefix + '.noIUPAC.' + sufix
        print prefix
        print args.input
        print temp
        print final
        convert_IUPAC  = '/rhome/cjinfeng/software/tools/fastatools/fastascripts/faunamb ' + args.input + ' > ' + temp
        convert_upcase = 'perl /rhome/cjinfeng/software/bin/fastaDeal.pl --reform upperize ' + temp + ' > ' + final
        os.system(convert_IUPAC)
        os.system(convert_upcase) 

if __name__ == '__main__':
    main()

