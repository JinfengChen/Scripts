perl /share/raid12/chenjinfeng/tools/bin/getidseq.pl -l PF03101_BR_63 -f /share/raid12/chenjinfeng/seqlib/BGI_analysis_data/bradi/Bradi_1.0.pep.fa -o PF03101_BR_63.fa
perl /share/raid12/chenjinfeng/tools/bin/getidseq.pl -l PF03101_OB_4 -f /share/raid12/chenjinfeng/seqlib/BGI_analysis_data/brachyantha/Oryza_brachyantha.genome.super_scaffold.v1.0.glean.gff.150_filter.2K.pep -o PF03101_OB_4.fa
perl /share/raid12/chenjinfeng/tools/bin/getidseq.pl -l PF03101_OS_47 -f /share/raid12/chenjinfeng/seqlib/BGI_analysis_data/rice_data/representative_orf.fa.nr.pep  -o PF03101_OS_47.fa


perl /share/raid12/chenjinfeng/tools/bin/getidseq.pl -l PF03101_OS_47 -f /share/raid12/chenjinfeng/seqlib/BGI_analysis_data/rice_data/representative_orf_nuc.fa.nr.cds -o PF03101_OS_47_cds.fa

blastall -i PF03101_OS_47_cds.fa -d /share/raid12/chenjinfeng/FFgenome/evolution/Bigfamily/input/super_scaffold.fa -p blastn -o PF03101_OS_47_cds2FF.blastm8 -m 8 -e 1e-5 > log 2> log2 &

