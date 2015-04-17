##bwa, short reads using aln and samse/sampe
##bwa, long reads using bwasw ( >200bp)
echo "fasta2fastq"
perl fasta2fastq.pl brachyantha_CenChIP.fna brachyantha_CenChIP.qual > brachyantha_CenChIP.fastq &
echo "index genome"
bwa index -a is OBa.chr.fa > log 2> log2 &
echo "map 454 reads to genome"
bwa bwasw OBa.chr.fa brachyantha_CenChIP.fastq > Cent4542OBa.sam &
echo "sam2bed"
perl /home/jfchen/FFproject/tools/bin/sam2bed.pl Cent4542OBa.sam Cent4542OBa.bed > log 2> log1 &
echo "sort"
msort -k1,n2 Cent4542OBa.bed > Cent4542OBa.sort.bed
perl bed2windowsbar.pl -win 200 -bar cent454.bar -bed Cent4542OBa.10.bed --sorting


echo "ffcentf to genome"
blastall -p blastn -i ffcent.txt -d OBa.chr.fa -o ffcent2OBa.allblast -e 1e-5 -m 8 > log &
echo "blast2bed"
perl blast2bed.pl --blastm8 ffcent2OBa.blast

