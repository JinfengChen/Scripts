find /home/jfchen/FFproject/FFgenome/05.epigenome/heterchromatin/HMMSeg/input/examples/ | grep Rna > rna.list
java -jar HMMSeg.jar --num-states 2 --input-bed --smooth 64000 --num-starts 10 --log log h3k27me3.list > log1 2> log2 &

perl bed2wave4methylation_2.pl --window 1000 --step 1000 --chrlen ../input/IRGSP.chrlen --bed ../input/rice_methylation_bed --bar ../input/rice_methylation_wave > log 2> log2 &
perl fillgap.pl --window 10000 -bed ../input/rice_methylation_wave/
find /home/jfchen/FFproject/FFgenome/05.epigenome/heterchromatin/HMMSeg/input/rice_methylation_wave | grep "chr02.CHG.wave.fill.bed" > chr02.chg10000.list
java -jar HMMSeg.jar --num-states 2 --input-bed --smooth 64000 --num-starts 10 --log log chr02.chg10000.list > log1 2> log2 &


perl bed2wave4methylation_2.pl --window 50000 --step 50000 --chrlen ../input/IRGSP.chrlen --bed ../input/rice_methylation_bed_test --bar ../input/rice_methylation_wave > log 2> log2 &
perl fillgap.pl --window 50000 -bed ../input/rice_methylation_wave/
find /home/jfchen/FFproject/FFgenome/05.epigenome/heterchromatin/HMMSeg/input/rice_methylation_wave | grep "chr02.CHG.wave.fill.bed" > chr02.chg50000.list
find /home/jfchen/FFproject/FFgenome/05.epigenome/heterchromatin/HMMSeg/input/rice_methylation_wave | grep "chr02.CG.wave.fill.bed" > chr02.cg50000.list
find /home/jfchen/FFproject/FFgenome/05.epigenome/heterchromatin/HMMSeg/input/rice_methylation_wave | grep "chr02.CHH.wave.fill.bed" > chr02.chh50000.list
java -jar HMMSeg.jar --num-states 2 --input-bed --smooth 200000 --num-starts 10 --log log chr02.chg10000.list > log1 2> log2 &
perl wig2bar.pl -wig ../input/rice_methylation_wave

perl fillgap.pl --window 50000 -bed ../input/rice_methylation_wave/


echo "OBa methylation"
perl bed2wave4methylation_2.pl -chrlen ../input/OBa.chrlen -bar ../input/OBa_methylation_wave -bed ../input/OBa_methylation_bed > log 2> log2 &
perl fillgap.pl --chrlen ../input/OBa.chrlen --bed ../input/OBa_methylation_wave --project OBa > log 2> log2 &

echo "OBa TE"
perl gff2wave4TE.pl --chrlen ../input/OBa.chrlen --bar ../input/OBa_feature_wave --gff ../input/OBa.all.fa.RepeatMasker.out.gff.chr.type > log 2> log2 &
perl smoothBED.pl --chrlen ../input/OBa.chrlen --bed ../input/OBa_feature_wave --project OBa > log 2> log2 &

echo "OBa manual TE"
perl gff2wave4TE.pl --chrlen ../input/OBa.chrlen --bar ../input/OBa_manual_feature_wave --gff ../input/OBa.all.manual.TE.gff.chr.type > log 2> log2 &
perl smoothBED.pl --chrlen ../input/OBa.chrlen --bed ../input/OBa_manual_feature_wave --project OBa > log 2> log2 &

echo "OBa absolute methylation"
perl bed2wave4methylation_2.pl --chrlen ../input/OBa.chrlen --bar ../input/OBa_absolute_methylation_wave --bed ../input/OBa_methylation_bed --absolute > log 2> log2 &
perl fillgap.pl --chrlen ../input/OBa.chrlen --bed ../input/OBa_absolute_methylation_wave --project OBa > log 2> log2 &


echo "IRGSP TE"
perl gff2wave4TE.pl --chrlen ../input/IRGSP.chrlen --bar ../input/IRGSP_feature_wave --gff ../input/IRGSP.build5.RepeatMasker.out.gff.chr.type > log 2> log2 &
perl smoothBED.pl --chrlen ../input/IRGSP.chrlen --bed ../input/IRGSP_feature_wave --project IRGSP > log 2> log2 &

echo "rice methylation"
perl bed2wave4methylation_2.pl --chrlen ../input/IRGSP.chrlen --bed ../input/rice_methylation_bed --bar ../input/rice_methylation_wave > log 2> log2 &
perl fillgap.pl --chrlen ../input/IRGSP.chrlen --bed ../input/rice_methylation_wave --project rice > log 2> log2 &

echo "rice gene"
perl gff2wave4GENE.pl --chrlen ../input/IRGSP.chrlen --bar ../input/rice_gene_wave --gff ../input/RAP3.gff3.nr.gff.chr > log 2> log2 &
perl smoothBED.pl --chrlen ../input/IRGSP.chrlen --bed ../input/rice_gene_wave --project rice

echo "OBa gene"
perl gff2wave4GENE.pl --chrlen ../input/OBa.chrlen --bar ../input/OBa_gene_wave --gff ../input/OBa.all.gff.chr > log 2> log2 &
perl smoothBED.pl --chrlen ../input/OBa.chrlen --bed ../input/OBa_gene_wave --project OBa

echo "TIGR6 H3K9me2"
perl bed2wave4chipseq.pl --chrlen ../input/TIGR6.chrlen --bar ../input/TIGR6_H3K9me2_wave --bed ../input/H3K9me2.bed.chr > log 2> log2 &
perl fillgap4chipseq.pl --chrlen ../input/TIGR6.chrlen --bed ../input/TIGR6_H3K9me2_wave --project TIGR6 > log 2> log2 &

