perl sumCm.pl -ref /home/jfchen/epigenome/peak/input/tigr6.1/all.fa -1 ../SRR042638_1.fastq -2 ../SRR042638_2.fastq -project BSmethyl > log &
perl step2_split2chr.pl -ref /home/jfchen/epigenome/peak/input/tigr6.1/all.fa -1 ../SRR042638_1.fastq -2 ../SRR042638_2.fastq -project BSmethyl > log 2>log2 &
