bowtie -v 2 -S ./bsgenome/Tgenome.fa -1 ./bsreads/SRR042638_1.fastq.T.fastq -2 ./bsreads/SRR042638_2.fastq.A.fastq SRR042638_1.fastq.TandSRR042638_2.fastq.A.T.sam > 1T2AtoT.log &
bowtie -v 2 -S ./bsgenome/Agenome.fa -1 ./bsreads/SRR042638_1.fastq.T.fastq -2 ./bsreads/SRR042638_2.fastq.A.fastq SRR042638_1.fastq.TandSRR042638_2.fastq.A.A.sam > 1T2AtoA.log &
