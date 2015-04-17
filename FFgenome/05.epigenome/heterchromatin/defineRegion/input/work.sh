perl ../bin/bed2windowsbar.pl -win 50 -bar CHGmethylation.chr04.bar.txt -bed ./rice_methylation_bed/chr04.CHG.bed > log &
perl /home/jfchen/FFproject/tools/bin/fastaDeal.pl --attr id:len /home/jfchen/FFproject/seqlib/BGI_analysis_data/rice_data/IRGSP.build5 > IRGSP.chrlen &
perl ../bin/bed2windowsbar4methylation.pl -win 50000 -bar CHGmethylation.chr04.bar.txt -bed ./rice_methylation_bed/chr04.CHG.bed > log &

perl ../bin/bed2windowsbar4methylation_2.pl -chrlen IRGSP.chrlen -bar rice_methylation_bar -bed rice_methylation_bed
perl ../bin/gff2windowsbar4gene.pl -chrlen IRGSP.chrlen -bar rice_gene_bar -gff RAP3.gff3.nr.gff.chr > log 2> log2 &
perl ../bin/gff2windowsbar4TE.pl -chrlen IRGSP.chrlen -bar rice_repeat_bar -gff IRGSP.build5.RepeatMasker.out.gff.chr > log 2> log2 &
perl ../bin/bed2windowsbar4chipseq.pl -chrlen IRGSP.chrlen -bar rice_chipseq_bar -type H3K4 -bed OS.H3K4.bed.chr > log 2> log2 &

perl ../bin/bed2windowsbar4methylation_2.pl --chrlen OBa.chrlen --bar OBa_methylation_bar -bed OBa_methylation_bed


echo "cent0 for FF"
perl /home/jfchen/FFproject/tools/bin/fastaDeal.pl --attr id:len /home/jfchen/FFproject/FFgenome/02.annotation/04.centromere/input/OBa.chr.fa > OBa.add100.chrlen
perl /home/jfchen/FFproject/tools/bin/splitGFF.pl Cent4542OBa.sort.bed
perl ../bin/bed2windowsbar4chipseq.pl --chrlen OBa.add100.chrlen --bar Cent4542OBa.sort.bed.bar --bed Cent4542OBa.sort.bed.chr --type CentO > log 2> log2 &

echo "centff for FF"
perl /home/jfchen/FFproject/tools/bin/splitGFF.pl ffcent2OBa.bed  
perl ../bin/bed2windowsbar4chipseq.pl --chrlen OBa.add100.chrlen --bar ffcent2OBa.bed.bar --bed ffcent2OBa.bed.chr --type CentF


