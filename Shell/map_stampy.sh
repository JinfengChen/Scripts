#!/bin/sh
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00

cd $PBS_O_WORKDIR

#ref=/rhome/cjinfeng/HEG4_cjinfeng/MappingReads/input/MSU_r7.fa
ref=/rhome/cjinfeng/HEG4_cjinfeng/Variations/mPing/mPing_reads/input/mping.fa
read1=HEG4_2.3.mPing.20X_p1.sorted.fq
read2=HEG4_2.3.mPing.20X_p2.sorted.fq
prefix=mPing.stampy.mping
#index=rice_MSU7
index=rice_mPing

echo "Mapping"
/opt/stampy/1.0.21-py2.7/stampy.py --species=rice --assembly=MSU7 -G $index $ref
/opt/stampy/1.0.21-py2.7/stampy.py -g $index -H $index
/opt/stampy/1.0.21-py2.7/stampy.py -g $index -h $index -M $read1 $read2 > $prefix.sam

echo "Convert Bam"
/usr/local/bin/samtools view -bS -o $prefix.raw.bam $prefix.sam
/usr/local/bin/samtools sort $prefix.raw.bam $prefix.sort
echo "RM dup"
java -jar /opt/picard/1.81/MarkDuplicates.jar ASSUME_SORTED=TRUE REMOVE_DUPLICATES=TRUE VALIDATION_STRINGENCY=LENIENT INPUT=$prefix.sort.bam OUTPUT=$prefix.bam METRICS_FILE=$prefix.dupli
echo "Clean"
#rm $prefix.sam $prefix.raw.bam $prefix.sort.bam
echo "Done

