#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=100:00:00

cd $PBS_O_WORKDIR

/usr/local/bin/samtools view out.bam -b -q 30 > out.unique.bam

echo "Done"

