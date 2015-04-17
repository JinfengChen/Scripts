 perl ../bin/protein_map_genome.pl  -cpu 100 -verbose -lines 100 -step 1234 ../input/rice.frags.fa.bgf.pep ../input/rice.frags.fa  

perl ../bin/protein_map_genome.pl -cpu 3 -verbose -step 1234 ../input/brachyantha_Fbox_Final.pep ../input/Gramene.chr.fa > log 2> log2 &

perl ../bin/protein_map_genome.pl -cpu 3 -verbose -step 1234 ../input/rice.frags.fa.bgf.pep ../input/rice.frags.fa > log 2> log2 &


