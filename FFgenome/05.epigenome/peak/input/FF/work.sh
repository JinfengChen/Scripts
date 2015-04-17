perl ../../bin/GFF2refFLAT.pl -g OBa.all.gff -r FFrefFlat.txt > log &
msort -k 3,n5 FFrefFlat.txt > FFrefFlat.txt.sort

cat /home/jfchen/FFproject/FFgenome/06.draw/drawFeature/data/gff/OBa.all.gff.chr/chr* > OBa.chr.gff
perl ../../bin/GFF2refFLAT.pl -g OBa.chr.gff -r OBarefFlat.txt
msort -k 3,n5 OBarefFlat.txt > OBarefFlat.txt.sort

