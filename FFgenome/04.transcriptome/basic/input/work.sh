echo "ff RPKM"
perl ../bin/GFF2refFLAT.pl -g Gramene.chr.gff -r Gramene.chr.refGene.txt
python ../bin/rpkm4gene.py -i Gramene.RNAseq.bowtie -a Gramene.chr.refGene.txt -o Gramene.gene.RPKM -bowtie -genePred -readcount
echo "tigr rice"
perl ../bin/GFF2refFLAT.pl -g tigr.all.final.gff -r tigr.all.final.refGene.txt
python ../bin/rpkm4gene.py -i tigr6.RNAseq.bowtie -a tigr.all.final.refGene.txt -o tigr6.gene.RPKM -bowtie -genePred -readcount



