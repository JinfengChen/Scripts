formatdb -i Gramene.chr.fa -p F
blastall -p blastn -i ./BAC/BAC.fasta -d Gramene.chr.fa -o BAC2Gramenev1.4.blastm8 -e 1e-5 -m 8 > log 2> log2 &
blastall -p blastn -i ./BAC/OB_B_PM_v6.txt -d Gramene.chr.fa -o Cent2Gramenev1.4.blastm8 -e 1e-5 -m 8 > log 2> log2 &
blastall -p blastn -i ./BAC/OB_B_PM_v6.txt -d Oryza_brachyantha.genome.super_scaffold.v1.0.fa -o Cent2Scaffold.blastm8 -e 1e-5 -m 8 > log 2> log2 &

blastall -p blastn -i ./chr1112/chr1112.fasta -d Gramene.chr.fa -o chr11122Gramenev1.4.blastm8 -e 1e-5 -m 8 > log 2> log2 &

blastall -p blastn -i ./chr03s/brachyantha.fa -d Gramene.chr.fa -o chr03s2Gramenev1.4.blastm8 -e 1e-5 -m 8 > log 2> log2 &

echo "set -a 10 to use 10 parrall running and set -n T to use MAGEblast, this is faster than defaut"
blastall -p blastn -a 10 -n T -i ./chr03s/brachyantha.fa -d Gramene.chr.fa -o chr03s2Gramenev1.4.blastm8 -e 1e-5 -m 8 > log 2> log2 &

blastall -p blastn -i ./chr03s/brachyantha.chr03s.fa -d Gramene.chr.fa -o chr03s2Gramenev1.4.blastm8 -e 1e-5 -m 8 > l g 2> log2 &





