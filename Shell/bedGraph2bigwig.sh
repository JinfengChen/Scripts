#!/bin/bash
 
# check commands: slopBed, bedGraphToBigWig and bedClip
 
which slopBed &>/dev/null || which bedtools &>/dev/null || { echo "slopBed/bedtools not found! Download bedTools: <http://code.google.com/p/bedtools/>"; exit 1; }
which bedGraphToBigWig &>/dev/null || { echo "bedGraphToBigWig not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }
which bedClip &>/dev/null || { echo "bedClip not found! Download: <http://hgdownload.cse.ucsc.edu/admin/exe/>"; exit 1; }
 
# end of checking
 
if [ $# -lt 2 ];then
    echo "Need 2 parameters! <bedgraph> <chrom info>"
    exit
fi
 
F=$1
G=$2
 
slopBed -i ${F} -g ${G} -b 0 | /rhome/cjinfeng/software/tools/UCSC_bin/bedClip stdin ${G} ${F}.clip
 
bedGraphToBigWig ${F}.clip ${G} ${F/bdg/bw}
 
rm -f ${F}.clip

