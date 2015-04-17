cat all_vs_all.blast.m8.solar.forHC.hcluster.stat | cut -f 3,4,5,6,7 > hcluster.matrix


echo "compare the cluster result from 159run and BGI"
cut -f3,4,5,6,7,8 ./TIGR6_159/all_vs_all.blast.m8.solar.forHC.hcluster.stat > ./TIGR6_159/tigr6_159.matrix
perl sumTable.pl --table ./TIGR6_159/tigr6_159.matrix > ./TIGR6_159/tigr6_159.matrix.table

cut -f3,4,5,6,7,8 ./BGI/all_vs_all.blast.m8.solar.forHC.hcluster.stat > BGI/bgi_IRGSP.matrix
perl sumTable.pl --table ./BGI/bgi_IRGSP.matrix > ./BGI/bgi_IRGSP.matrix.table

echo "compare the cluster result from gramene and BGI"
cut -f3,4,5 ./TIGR6_GRAMENE/all_vs_all.blast.m8.solar.forHC.hcluster.stat > ./TIGR6_GRAMENE/gramene.matrix
perl sumTable.pl --table ./TIGR6_GRAMENE/gramene.matrix > ./TIGR6_GRAMENE/gramene.matrix.table

echo "gramene_v1_final6_best2"
cp /home/jfchen/FFproject/FFgenome/03.evolution/ortholog_paralog_family/output/all_vs_all.blast.m8.solar.forHC.hcluster.stat TIGR6_GRAMENE_v1_final6_best2/
cut -f3,4,5 TIGR6_GRAMENE_v1_final6_best2/all_vs_all.blast.m8.solar.forHC.hcluster.stat > TIGR6_GRAMENE_v1_final6_best2/gramene.matrix
perl sumTable.pl --table TIGR6_GRAMENE_v1_final6_best2/gramene.matrix > TIGR6_GRAMENE_v1_final6_best2/gramene.matrix.table

echo "summary gene lose and family size changes"
perl sumstat.pl --table TIGR6_OBa_SB/all_vs_all.blast.m8.solar.forHC.hcluster.stat > summary
perl getidipr.pl -l OBRACH_lose.gene -i /home/jfchen/FFproject/seqlib/BGI_analysis_data/tigr_rice/tigr.all.final.iprscan -o OBRACH_lose.gene.iprscan
perl top50iprscan.pl --iprscan OBRACH_lose.gene.iprscan > OBRACH_lose.gene.iprscan.sum

perl getidipr.pl -l TIGR6_lose.gene -i /home/jfchen/FFproject/seqlib/BGI_analysis_data/brachyantha/superscaffold/Oryza_brachyantha.genome.super_scaffold.v1.0.glean.gff.150_filter.2K.pep.iprscan -o TIGR6_lose.gene.iprscan
perl top50iprscan.pl --iprscan TIGR6_lose.gene.iprscan > IGR6_lose.gene.iprscan.sum

perl NonsyntenyFunction.pl --gff /home/jfchen/FFproject/seqlib/BGI_analysis_data/tigr_rice/tigr.all.final.gff --id /ho
me/jfchen/FFproject/FFgenome/03.evolution/genefamily_hcluster/sumType/bin/TIGR6_3way/OrthOrNot/OBRACH_TIGR6/TIGR6.nono
rtholog.table > log

perl NonsyntenyFunction.pl --gff /home/jfchen/FFproject/seqlib/BGI_analysis_data/tigr_rice/tigr.all.final.gff --id /ho
me/jfchen/FFproject/FFgenome/03.evolution/genefamily_hcluster/sumType/bin/TIGR6_3way/OrthOrNot/OBRACH_TIGR6/TIGR6.nono
rtholog.table --iprscan /home/jfchen/FFproject/seqlib/BGI_analysis_data/tigr_rice/tigr.all.final.iprscan > log


perl NonsyntenyFunction.pl --gff /home/jfchen/FFproject/seqlib/BGI_analysis_data/tigr_rice/tigr.all.final.gff --id /ho
me/jfchen/FFproject/FFgenome/03.evolution/genefamily_hcluster/sumType/bin/TIGR6_3way/ClusterGeneType/TIGR6.unique.tabl
e --iprscan /home/jfchen/FFproject/seqlib/BGI_analysis_data/tigr_rice/tigr.all.final.iprscan > log1

perl NonsyntenyFunction.pl --gff /home/jfchen/FFproject/seqlib/BGI_analysis_data/tigr_rice/tigr.all.final.gff --id /ho
me/jfchen/FFproject/FFgenome/03.evolution/genefamily_hcluster/sumType/bin/TIGR6_3way/ClusterGeneType/TIGR6.noncluster.
table --iprscan /home/jfchen/FFproject/seqlib/BGI_analysis_data/tigr_rice/tigr.all.final.iprscan > log2

