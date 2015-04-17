perl /home/jfchen/FFproject/tools/bin/getidseq.pl --l anchorgene.list --f /home/jfchen/FFproject/seqlib/BGI_analysis_data/tigr_rice/tigr.all.pep.final.fa -o anchor.pep.fa

blastall -p tblastn -i anchor.pep.fa -d /home/jfchen/FFproject/seqlib/BGI_analysis_data/Sorghum_bicolor/Sb.genome.fa -o anchor2sorghum.blastm8 -m 8 -e 1e-5 > log 2> log2 &
blastall -p tblastn -i anchor.pep.fa -d /home/jfchen/FFproject/seqlib/BGI_analysis_data/bradi/Brachypodium_distachyon_Bd21.main_genome.scaffolds.fasta -o anchor2bradi.blastm8 -m 8 -e 1e-5 > log 2> log2 &


perl GetSubGFF.pl --qryhead OS_chr12_24586000_25020000 --qrygff3 /home/jfchen/FFproject/seqlib/BGI_analysis_data/tigr_rice/tigr.all.final.gff --qryfasta /home/jfchen/FFproject/seqlib/BGI_analysis_data/tigr_rice/all.con --project OS 

perl GetSubGFF.pl --qryhead OB_chr12_13158000_13450000 --qrygff3 /home/jfchen/FFproject/seqlib/BGI_analysis_data/brachyantha/gramene/input/Gramenev1.4/Gramene.chr.gff --qryfasta /home/jfchen/FFproject/seqlib/BGI_analysis_data/brachyantha/gramene/input/Gramenev1.4/Gramene.chr.fa --project OB &

perl GetSubGFF.pl --qryhead BD_Bd4_2044000_2474000 --qrygff3 /home/jfchen/FFproject/seqlib/BGI_analysis_data/bradi/Bradi_1.0.final.gff --qryfasta /home/jfchen/FFproject/seqlib/BGI_analysis_data/bradi/Brachypodium_distachyon_Bd21.main_genome.scaffolds.fasta --project BD &

perl GetSubGFF.pl --qryhead SB_chromosome8_50771000_51195000 --qrygff3 /home/jfchen/FFproject/seqlib/BGI_analysis_data/Sorghum_bicolor/Sb.final.gff --qryfasta /home/jfchen/FFproject/seqlib/BGI_analysis_data/Sorghum_bicolor/Sb.genome.fa --project SB &



perl runblast2seq.pl
perl run2act.pl
