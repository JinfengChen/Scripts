bowtie-build ../input/LargeVerifySeq.fas LargeVerifySeq >log 2>log2
tophat -r 50 LargeVerifySeq ../fq/Root_R8_L1_1.fq ../fq/Root_R8_L1_2.fq > log 2> log2 &

bowtie-build ../input/OBa.all.fa OBa.all > log 2> log2 &
perl ../bin/runtophat.pl -ref ./OBa.all -length ./reflen.txt -1 ../fq/Root_1 -2 ../fq/Root_2 -p OBaTest

perl ../bin/runtophat.pl -ref ./OBa.all -length ./reflen.txt -1 ../fq/Root_R8_L1_1.fq -2 ../fq/Root_R8_L1_2.fq -p OBaRoot > log 2> log2 &

sed s/CDS/exon/g /share/raid12/chenjinfeng/FFgenome/assmbly/superscaf/bin/OBa.all.gff > OBa.all.exon.gff

perl ../bin/runtophat.pl -ref ./OBa.all -length ./reflen.txt -gff ./OBa.all.exon.gff -1 ../fq/Root_R8_L1_1.fq -2 ../fq/Root_R8_L1_2.fq -p OBaRoot_GFF > log 2> log2  &


perl ../bin/runtophat.pl -ref ./30 -length ./30.len -gff ./30.exon.gff -1 ../fq/Root_R8_L1_1.fq -2 ../fq/Root_R8_L1_2.fq -p 30Root_GFF > log 2> log2 &

perl PrepareGFF.pl -gff /share/raid12/chenjinfeng/FFgenome/assmbly/superscaf/bin/OBa.all.gff

perl ../bin/runtophat.pl -ref ./OBa.all -length ./reflen.txt -gff ./OBa.all.tophat.gff -1 ../fq/Root_1.fq -2 ../fq/Root_2.fq -p OBaRoot_GFF > OBaRoot_GFF.log 2> OBaRoot_GFF.log2  &

perl ../bin/runtophat.pl -ref ./IRGSP.build5 -length ./IRGSPlen.txt -gff ./RAP3.gff3.nr.tophat.gff -read ../fq/shoot_japonica.fq -p IRGSP_Shoot_GFF > log 2> log2 &

perl ../bin/runtophat.pl -ref ./OBa.all -length ./reflen.txt -gff ./OBa.all.tophat.gff -1 ../fq/shoot_2_L1_1.fq -2 ../fq/shoot_2_L1_2.fq -p OBaShoot_GFF > OBaShoot_GFF.log 2> OBaShoot_GFF.log2  &


perl ../bin/runtophat.pl -ref ./OBa.all -length ./reflen.txt -gff ./OBa.all.tophat.gff -1 /share/raid12/chenjinfeng/FFgenome/transcriptome/fq/shoot_2_L1_1.fq,/share/raid12/chenjinfeng/FFgenome/transcriptome/fq/Root_R8_L1_1.fq -2 /share/raid12/chenjinfeng/FFgenome/transcriptome/fq/shoot_2_L1_2.fq,/share/raid12/chenjinfeng/FFgenome/transcriptome/fq/Root_R8_L1_2.fq -p OBa_GFF_alljunction > OBa_GFF.log 2> OBa_GFF.log2 &


perl RNAseqDistribution.pl -bed OBaShoot_tophat/accepted_hits.bed -gene OBa.all.gff -repeat OBa.all.fa.RepeatMasker.out.gff -project testsum > log 2> log2 &

perl RNAseqDistribution.pl -bed testsum.bed -gene testsum.gene.gff -repeat testsum.repeat.gff -project testsum > log 2> log2 &

