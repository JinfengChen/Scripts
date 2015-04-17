#!/opt/Python/2.7.3/bin/python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import argparse

def usage():
    message='Python ChrName.py --input INPUT --output OUTPUT'
    print message

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output')
    parser.add_argument('-v', dest='verbose', action='store_true')
    args = parser.parse_args()
    try:
        len(args.output) > 1
    except:
        usage()
        sys.exit(2)

    ofile = open(args.output,"w")
    for record in SeqIO.parse(args.input,"fasta"):
        print record.id
        id = re.match('\D*(\d+)',record.id)
        chr = "chr" + id.group(1).zfill(2)
        newrecord = SeqRecord(record.seq,id=chr,description="")
        SeqIO.write(newrecord,ofile,"fasta")
    ofile.close()

