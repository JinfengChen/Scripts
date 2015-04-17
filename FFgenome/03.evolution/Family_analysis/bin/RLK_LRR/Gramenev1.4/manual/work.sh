perl RLKLRR.pl --ortholog ortholog.txt --tandem tandem.txt
perl Typelist.pl --dir ./type/ > Gramenev1.4.RLK_LRR.ID
perl checkID2.pl --id1 manaul.RLKLRR.ID --id2 ../Gramenev1.4.RLK_LRR.ID
perl ../../../getidpos.pl --list miss.list --gff /home/jfchen/FFproject/seqlib/BGI_analysis_data/brachyantha/gramene/input/Gramenev1.4/Gramene.chr.gff --output miss.list.inf


perl RLKLRR.pl --ortholog RKL_LRR_final_ortholog.txt --tandem RKL_LRR_final_tandem.txt --transpose1 RKL_LRR_final_transposed_rice.txt --transpose2 RKL_LRR_final_transposed_FF.txt

perl RLKLRR.pl --ortholog RKL_LRR_final_ortholog.txt --tandem RKL_LRR_final_tandem.txt --transpose1 RKL_LRR_final_transposed_rice.txt --transpose2 RKL_LRR_final_transposed_FF.txt > final.summary

echo "pseudogene in FF may have functional gene in rice, these gene were added into final list manually"

perl gettypeid.pl --id manaul.rice.RLKLRR.ID --type LRR-XII > LRR_XII.rice.id
perl gettypeid.pl --id manaul.OBa.RLKLRR.ID --type LRR-XII > LRR_XII.OBa.id
perl ../../../getidseq.pl --list LRR_XII.rice.id --fasta ../../../../input/tigr.all.pep.final.fa --output LRR_XII.rice.fa
perl ../../../getidseq.pl --list LRR_XII.OBa.id --fasta ../../../../input/Gramene.pep.fa --output LRR_XII.OBa.fa
cat LRR_XII.rice.fa LRR_XII.OBa.fa > LRR_XII.fa
perl /home/jfchen/FFproject/FFgenome/03.evolution/Family_analysis/bin/DrawTree.pl --protein LRR_XII.fa --align --nj > log 2> log2 &


