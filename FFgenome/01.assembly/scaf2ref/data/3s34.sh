/share/raid1/genome/bin/formatdb -i chr34.con -p F
/share/raid1/genome/bin/blastall -p blastn -i brachyantha_3S_scaffolds.fna -d chr34.con -o 3s2chr34blastm8 -e 1e-5 -m 8 > log &
