perl ../../bin/GFF2refFLAT.pl -g RAP3.gff3.nr.gff -r IRGSPrefFlat.txt > log &
msort -k 3,n5 IRGSPrefFlat.txt > IRGSPrefFlat.txt.sort
