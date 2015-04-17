echo "FF"
bowtie-build Gramene.chr.fa Gramene.chr.fa
bowtie Gramene.chr.fa -1 ../fq/Root_R8_L1_1.fq,../fq/shoot_2_L1_1.fq -2 ../fq/Root_R8_L1_2.fq,../fq/shoot_2_L1_2.fq  Gramene.RNAseq.bowtie > log 2> log2 &

echo "tigr rice"
bowtie-build all.con all.con > log 2> log2 &
bowtie all.con ../fq/shoot_japonica.fq tigr6.RNAseq.bowtie > log 2> log2 &


