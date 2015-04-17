#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=2g
#PBS -l walltime=100:00:00

#cd $PBS_O_WORKDIR

module load bamtools
module load samtools

bam=/rhome/cjinfeng/BigData/00.RD/MappingReads/bin/A119_A123allpathv1_BWA/Chromosome/A119.HEG4allpathv1_BWA.ALL.bam
#bamtools random -in $bam -out A119.HEG4allpathv1_BWA.ALL.random.bam -n 100000000
#samtools sort A119.HEG4allpathv1_BWA.ALL.random.bam A119.HEG4allpathv1_BWA.ALL.random.sort

java -Xmx2000m -jar DownsampleSam.jar I=$bam O=A119.HEG4allpathv1_BWA.ALL.random.bam P=0.2
#samtools sort A119.HEG4allpathv1_BWA.ALL.random.bam A119.HEG4allpathv1_BWA.ALL.random.sort


echo "Done"
