perl /home/jfchen/159/FFproject/tools/bin/fastaDeal.pl -sub 7741500-8931000 brachyantha.chr03s.fa > brachyantha.chr03s.miss.fa

perl /home/jfchen/159/FFproject/tools/bin/getsubgff.pl --refseq chr03 --start 7741500 --end 8931000 --gff brachyantha.chr03s.gene.gff --output brachyantha.chr03s.miss.gene.gff

perl /home/jfchen/159/FFproject/tools/bin/getGene.pl brachyantha.chr03s.miss.gene.gff brachyantha.chr03s.miss.fa > brachyantha.chr03s.miss.cds

perl /home/jfchen/159/FFproject/tools/bin/cds2aa.pl brachyantha.chr03s.miss.cds > brachyantha.chr03s.miss.pep

blastall -p blastp -i brachyantha.chr03s.miss.pep -d /home/jfchen/159/FFproject/seqlib/BGI_analysis_data/tigr_rice/all.pep -o miss2tigrpep.blastm8 -m 8 -e 1e-5 > log 2> log2 &

blastall -a 10 -p blastp -i brachyantha.chr03s.miss.pep -d /home/jfchen/159/FFproject/seqlib/BGI_analysis_data/tigr_rice/all.pep -o miss2tigrpep.blastm8 -m 8 -e 1e-5 > log 2> log2 &

blastall -a 10 -p blastp -i brachyantha.chr03s.miss.pep -d /home/jfchen/159/FFproject/seqlib/BGI_analysis_data/tigr_rice/all.pep -o miss2tigrpep.blast -e 1e-5 > log 2> log2 &

perl /home/jfchen/159/FFproject/tools/bin/blast_parser_bgi.pl --tophit 1 miss2tigrpep.blast > miss2tigrpep.blast.table




