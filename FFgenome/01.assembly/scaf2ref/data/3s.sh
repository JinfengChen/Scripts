formatdb -i chr03.con -p F
blastall -p blastn -i brachyantha_3S_scaffolds.fna -d chr03.con -o 3s2chr03blastm8 -e 1e-5 -m 8 > log &
