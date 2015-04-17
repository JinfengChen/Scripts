#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l mem=10gb
#PBS -l walltime=10:00:00

cd $PBS_O_WORKDIR

for i in `ls *.group.bam | sed 's/@//'`;
do
   input+="INPUT=$i "
done

dirname=`pwd`
prefix=`basename $dirname`

echo $input
echo $prefix

#java -jar /opt/picard/1.81/MergeSamFiles.jar VALIDATION_STRINGENCY=LENIENT $input USE_THREADING=true OUTPUT=$prefix.merge.bam

echo "Done"
