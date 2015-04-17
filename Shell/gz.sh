#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=10:00:00

cd $PBS_O_WORKDIR


for i in `ls *.fastq | sed 's/@//'`
do
   echo $i
   if [ ! -e $i.gz ]; then
      gzip $i
   fi
done

