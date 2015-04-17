#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=10:00:00

cd $PBS_O_WORKDIR


for i in `ls *.bam | sed 's/@//'`
do
   echo $i
   prefix=basename $i .bam 
   if [ ! -e $i.bai ]; then
   /usr/local/bin/samtools sort $i $prefix.sort.bam
   /usr/local/bin/samtools index $prefix.sort.bam
   fi
done

