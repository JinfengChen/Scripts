bowtie -v 2 -S ./bsgenome/Tgenome.fa -1 ./bsreads/SRR042638_1.fastq.T.fastq -2 ./bsreads/SRR042638_2.fastq.A.fastq SRR042638_1.fastq.TandSRR042638_2.fastq.A.T.sam > 1T2AtoT.log &
bowtie -v 2 -S ./bsgenome/Agenome.fa -1 ./bsreads/SRR042638_1.fastq.T.fastq -2 ./bsreads/SRR042638_2.fastq.A.fastq SRR042638_1.fastq.TandSRR042638_2.fastq.A.A.sam > 1T2AtoA.log &
samtools view -bS -o 1T2A.A.bam SRR042638_1.fastq.TandSRR042638_2.fastq.A.A.sam
samtools view -bS -o 1T2A.T.bam SRR042638_1.fastq.TandSRR042638_2.fastq.A.T.sam
bamToBed -i 1T2A.A.bam > 1T2A.A.bed
bamToBed -i 1T2A.T.bam > 1T2A.T.bed
perl sumCm.pl -ref /home/jfchen/epigenome/peak/input/tigr6.1/all.fa -1 ../SRR042638_1.fastq -2 ../SRR042638_2.fastq -project BSmethyl > log
msort -k 1,n2 BSmethyl.plus.status > BSmethyl.plus.status.sorted
msort -k 1,n2 BSmethyl.minus.status > BSmethyl.minus.status.sorted
perl analyzeCm.pl -ref /home/jfchen/epigenome/peak/input/tigr6.1/all.fa -plus BSmethyl.plus.status.sorted -minus BSmethyl.minus.status.sorted > log &
perl analyzeCm.pl -ref /home/jfchen/epigenome/peak/input/tigr6.1/all.fa -plus BSmethyl.plus.status -minus BSmethyl.minus.status > log &
