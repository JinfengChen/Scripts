#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=10:00:00

cd $PBS_O_WORKDIR

for i in `ls *.bam | sed 's/@//'`
do
   echo $i
   (( j += 1 ))
   echo $j
   if [ ! -e $i.group.bam ]; then
   #/usr/local/bin/samtools reheader $i.header $i > $i.group.bam
   java -jar /opt/picard/1.81/AddOrReplaceReadGroups.jar VALIDATION_STRINGENCY=LENIENT I=$i O=$i.group.bam ID=$j LB=HEG4_2.3 PL=illumina PU=0$j SM=HEG4
   echo $j
   fi
done

