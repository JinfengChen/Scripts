echo "get gff from different cluster types, OBa"
perl getidgff.pl --l ../input/TIGR6_GRAMENE_v1_final8_3_last/ClusterGeneType/OBRACH.noncluster.table --g ../input/OBa.gramene.FINAL8_3_last.bestmodel.filter.gff --o FINAL8_3_last.noncluster.gff
perl getidgff.pl --l ../input/TIGR6_GRAMENE_v1_final8_3_last/ClusterGeneType/OBRACH.unique.table --g ../input/OBa.gramene.FINAL8_3_last.bestmodel.filter.gff --o FINAL8_3_last.unique.gff
perl getidgff.pl --l ../input/TIGR6_GRAMENE_v1_final8_3_last/ClusterGeneType/OBRACH.shared.table --g ../input/OBa.gramene.FINAL8_3_last.bestmodel.filter.gff --o FINAL8_3_last.share.gff

echo "get stat for each type"
perl gff_statisticV2.pl FINAL8_3_last.noncluster.gff final8_3.noncluster
perl gff_statisticV2.pl FINAL8_3_last.unique.gff final8_3.unique
perl gff_statisticV2.pl FINAL8_3_last.share.gff final8_3.share

echo "compare among types and draw"
perl CompareGeneStat.pl --project final8_3 > log 2> log2 &

echo "get gff from different cluster types, rice"
perl getidgff.pl --l ../input/TIGR6_GRAMENE_v1_final8_3_last/ClusterGeneType/TIGR6.noncluster.table --g ../input/tigr.all.final.gff --o tigr.noncluster.gff
perl getidgff.pl --l ../input/TIGR6_GRAMENE_v1_final8_3_last/ClusterGeneType/TIGR6.unique.table --g ../input/tigr.all.final.gff --o tigr.unique.gff
perl getidgff.pl --l ../input/TIGR6_GRAMENE_v1_final8_3_last/ClusterGeneType/TIGR6.shared.table --g ../input/tigr.all.final.gff --o tigr.share.gff

echo "get stat for each type"
perl gff_statisticV2.pl tigr.noncluster.gff tigr.noncluster
perl gff_statisticV2.pl tigr.unique.gff tigr.unique
perl gff_statisticV2.pl tigr.share.gff tigr.share

echo "compare among types and draw"
perl CompareGeneStat.pl --project tigr > log 2> log2 &

echo "compare between rice and FF"
perl Compare2SpeciesGeneStat.pl --project1 tigr --project2 final8_3

echo "count annotated gene"
perl CompareGeneIPR.pl --iprscan ../input/iprscan/tigr.all.final.iprscan --dir ../input/TIGR6_GRAMENE_v1_final8_3_last/ClusterGeneType --project TIGR6
perl CompareGeneIPR.pl --iprscan ../input/iprscan/OBa.gramene.FINAL8_3_last.bestmodel.filter.pep.iprscan --dir ../input/TIGR6_GRAMENE_v1_final8_3_last/ClusterGeneType --project OBRACH
perl CompareGeneIPR.pl --iprscan ../input/iprscan/Sb.final.pep.iprscan --dir ../input/TIGR6_GRAMENE_v1_final8_3_last/ClusterGeneType --project SBICO

echo "top 50 iprscan"
perl top50iprscan.pl --iprscan ../input/iprscan/tigr.all.final.iprscan --gff ../input/tigr.all.final.gff > tigr.top50 2> log2 &
perl top50iprscan.pl --iprscan ../input/iprscan/OBa.gramene.FINAL8_3_last.bestmodel.filter.pep.iprscan --gff ../input/OBa.gramene.FINAL8_3_last.bestmodel.filter.gff > OBa.top50 2> log2 &

echo "get subtype iprscan"
perl getidipr.pl -l ../input/TIGR6_GRAMENE_v1_final8_3_last/ClusterGeneType/TIGR6.noncluster.table -i ../input/iprscan/tigr.all.final.iprscan -o tigr.noncluster.iprscan

