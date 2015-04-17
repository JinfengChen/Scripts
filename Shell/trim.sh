#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=2gb
#PBS -l walltime=10:00:00

cd $PBS_O_WORKDIR

zcat ERR035777.fastq.gz | fastx_clipper -a TGGAATTCTCGGGTGCCAAGGC -l 12 -Q 33 -v -c -z -o ERR035777.3trim.fastq.gz

echo "Done"

