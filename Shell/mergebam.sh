#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=10:00:00

cd $PBS_O_WORKDIR


java -jar /opt/picard/1.81/MergeSamFiles.jar VALIDATION_STRINGENCY=LENIENT INPUT=FC52_7.MSU7_BWA.bam.group.bam INPUT=FC52_8.MSU7_BWA.bam.group.bam USE_THREADING=true OUTPUT=HEG4.merge.bam


