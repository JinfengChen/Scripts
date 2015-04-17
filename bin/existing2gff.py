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
python existing2gff.py --input EG4.mping.all_reference.txt

    '''
    print message

#strain  TE      existingTE_coor reads_align_2_start     reads_align_2_end
#EG4     mping   Chr1:24779771..24780200 0       0
def txt2gff(infile, outfile):
    #print infile, outfile
    ofile = open(outfile, 'w')
    count = 0
    r_pos = re.compile(r'(\w+):(\d+)\.\.(\d+)')
    with open (infile, 'r') as filehd:
        for line in filehd:
            line = line.rstrip()
            #print line
            if len(line) > 2 and not line.startswith(r'strain'):
                unit = re.split(r'\t',line)
                count += 1
                chro, start, end = ['', 0, 0]
                m = r_pos.search(unit[2])
                if m: 
                    chro  = m.groups(0)[0]
                    start = m.groups(0)[1]
                    end   = m.groups(0)[2]
                #r_count = re.sub(r'\D+', '', unit[7])
                #l_count = re.sub(r'\D+', '', unit[8])
                #r_supp  = re.sub(r'\D+', '', unit[10])
                #l_supp  = re.sub(r'\D+', '', unit[11])
                r_id    = 'repeat_%s_%s_%s' %(chro, start, end)
                print >> ofile, '%s\t%s\t%s\t%s\t%s\t.\t.\t.\tID=%s;' %(chro, 'Nipponbare', 'blat',start, end, r_id)
    ofile.close()


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
    
    txt2gff(args.input, 'Nipponbare.mPing.gff')


if __name__ == '__main__':
    main()

