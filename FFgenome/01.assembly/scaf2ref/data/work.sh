#formatdb -i chr04.txt -p F
blastall -p blastn -i chr04.scafSeq -d chr04.txt -o scaf2chr04blastm8 -e 1e-5 -m 8 > log &
