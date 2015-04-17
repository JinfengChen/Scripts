perl ../../../bin/bsmap.pl -ref /home/jfchen/epigenome/peak/input/tigr6.1/all.fa -1 SRR042638_1.fastq -2 SRR042638_2.fastq -project BSmethyl
perl ../../../bin/bsmap.pl -ref /home/jfchen/epigenome/peak/input/tigr6.1/all.fa -1 read1a.fastq,read1b.fastq -2 read2a.fastq,read2b.fastq -project BSmethylPE > log &
perl ../../../bin/bsmap.pl -ref /home/jfchen/epigenome/peak/input/tigr6.1/all.fa -read read1a.fastq,read1b.fastq -project BSmethylSE > log &
perl ../../../bin/step1_bsmap.pl -ref /home/jfchen/epigenome/peak/input/tigr6.1/all.fa -read SRR042639.fastq,SRR042640.fastq,SRR042641.fastq,SRR042642.fastq -project BSmethylSE > log &
