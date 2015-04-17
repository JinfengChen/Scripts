#!/usr/bin/env python

#
# Author:       Alex Reynolds
#
# Project:      Converts 1-based, closed [a,b] GFF3 input
#               into 0-based, half-open [a-1,b) six-column extended BED
#
# Version:      1.2.2k2
#
# Notes:        The GFF3 specification (http://www.sequenceontology.org/gff3.shtml) 
#               contains columns that do not map directly to common or UCSC BED columns.
#               Therefore, we add the following columns to preserve the ability to 
#               seamlessly convert back to GFF3 after performing operations with 
#               BEDOPS or other tools.
#
#               - The 'source' GFF column data maps to the 7th BED column
#               - The 'type' data maps to the 8th BED column
#               - The 'phase' data maps to the 9th BED column
#               - The 'attributes' data maps to the 10th BED column
#
#               We make the following assumptions about the GFF3 input data:
#
#               - The 'seqid' GFF column data maps to the chromosome label (1st BED column)
#               - The 'ID' attribute in the 'attributes' GFF column (if present) maps to 
#                 the element ID (4th BED column)
#               - The 'score' and 'strand' GFF columns (if present) are mapped to the
#                 5th and 6th BED columns, respectively
#
#               If we encounter zero-length insertion elements (which are defined
#               where the start and stop GFF column data values are equivalent), the 
#               start coordinate is decremented to convert to 0-based, half-open indexing, 
#               and a 'zeroLengthInsertion' attribute is added to the 'attributes' GFF 
#               column data.
#
# Example usage:
#
#  $ gff2bed.py < foo.gff | sort-bed -
#

import sys, os

def main(*args):
    requiredVersion = (2,5)
    checkInstallation(requiredVersion)

    for line in sys.stdin:
        chomped_line = line.rstrip(os.linesep)
        if chomped_line.startswith('##'):
            print chomped_line
        else:
            elems = chomped_line.split('\t')
            cols = dict()
            cols['seqid'] = elems[0]
            cols['source'] = elems[1]
            cols['type'] = elems[2]
            cols['start'] = int(elems[3])
            cols['end'] = int(elems[4])
            cols['score'] = elems[5]
            cols['strand'] = elems[6]
            cols['phase'] = elems[7]
            cols['attributes'] = elems[8]

            attrd = dict()
            attrs = map(lambda s: s.split('='), cols['attributes'].split(';'))
            for attr in attrs:
                attrd[attr[0]] = attr[1]

            cols['chr'] = cols['seqid']
            try:
                cols['id'] = attrd['ID']
            except KeyError:
                cols['id'] = '.'

            if cols['start'] == cols['end']:
                cols['start'] -= 1
                cols['attributes'] = ';'.join([cols['attributes'], "zeroLengthInsertion=True"])
            else:
                cols['start'] -= 1

            print '\t'.join([cols['chr'], 
                             str(cols['start']),
                             str(cols['end']),
                             cols['id'],
                             cols['score'],
                             cols['strand'],
                             cols['source'],
                             cols['type'],
                             cols['phase'],
                             cols['attributes']])
    return 0

def checkInstallation(rv):
    currentVersion = sys.version_info
    if currentVersion[0] == rv[0] and currentVersion[1] >= rv[1]:
        pass
    else:
        sys.stderr.write( "[%s] - Error: Your Python interpreter must be %d.%d or greater (within major version %d)\n" % (sys.argv[0], rv[0], rv[1], rv[0]) )
        sys.exit(-1)
    return 0

if __name__ == '__main__':
    sys.exit(main(*sys.argv))


