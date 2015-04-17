#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=10:00:00

cd $PBS_O_WORKDIR

BAM=HEG4_2.3.MSU7_BWA.bam
#BAM=out.bam
OUT=HEG4_2.3.MSU7_BWA
for chr in `seq 1 12`
do
    echo Chr$chr
    samtools view -bh $BAM Chr${chr} | samtools sort - $OUT.Chr${chr}
    samtools index $OUT.Chr${chr}.bam
done
