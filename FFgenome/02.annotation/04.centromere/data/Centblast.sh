blastall -p blastn -i CentO.txt -d /share/raid12/chenjinfeng/seqlib/rice/all.con -o CentO2riceblast -e 1e-5 -m 8 >log&
blastall -p blastn -i CentO.txt -d /share/raid12/chenjinfeng/seqlib/BGI_analysis_data/rice_data/IRGSP.build5 -o CentO2IRASPblast -e 1e-5 -m 8 > log &
blastall -p blastn -i ffcent.txt -d /share/raid12/chenjinfeng/FFgenome/assmbly/superscaf/bin/OBa.all.fa -o ffcent2OBa.allblast -e 1e-5 -m 8 > log &
blastall -p blastn -i ffcent.txt -d /home/jfchen/FFproject/seqlib/BGI_analysis_data/brachyantha/gramene/input/Gramenev1.4/Gramene.chr.fa -o ffcent2Gramenev1.4blast -e 1e-5 -m 8 > log 2> log2 &
blastall -p blastn -i ffchrcent.txt -d /home/jfchen/FFproject/seqlib/BGI_analysis_data/brachyantha/gramene/input/Gramenev1.4/Gramene.chr.fa -o ffchrcent2Gramenev1.4blast -e 1e-5 -m 8 > log 2> log2 &


