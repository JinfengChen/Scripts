perl gtf2gff.pl < ../input/transcripts.gtf --gff3 --out transcripts.gff

perl getcuflink.pl ../input/transcripts.gff ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa > log 2> log2 &

echo "get first 100 sequence and blast with soapdenovo assembl to check"
perl /home/jfchen/FFproject/tools/bin/fastaDeal.pl --sample 1-100 cuflink.tran.fa > 100.fa
blastall -i 100.fa -d ../input/rice.100scaftig -p blastn -e 1e-5 -o 100blast > log 2> log2 &
