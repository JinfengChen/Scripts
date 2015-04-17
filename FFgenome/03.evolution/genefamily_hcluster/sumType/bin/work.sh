perl ClusterGeneType.pl --cluster ../input/all_vs_all.blast.m8.solar.forHC.hcluster --blastm8 ../input/all_vs_all.blast.m8 --cds ../input/all.cds > log 2> log2 &
./sumdistance.sh

echo "OBa rice"
perl OrthorNot.pl --typedir ./TIGR6/ClusterGeneType --distance ./TIGR6/Distance/OBRACH_TIGR6.distance.txt > log 2> log2 &
perl SyntenyOrNot.pl --syntid ./TIGR6/ob2os_TIGR6/ob_by_os.dist20.synteny_list --orthid ./TIGR6/OrthOrNot/OBRACH_TIGR6/TIGR6.ortholog.table --outdir ./TIGR6/SytnteyOrNot --project TIGR6
perl SyntenyOrNot.pl --syntid ./TIGR6/ob2os_TIGR6/os_by_ob.dist20.synteny_list --orthid ./TIGR6/OrthOrNot/OBRACH_TIGR6/OBRACH.ortholog.table --outdir ./TIGR6/SytnteyOrNot --project OBRACH

echo "OBa sorghum"
perl OrthorNot.pl --typedir ./TIGR6/ClusterGeneType --distance ./TIGR6/Distance/OBRACH_SBICO.distance.txt > log 2> log2 &
perl SyntenyOrNot.pl --syntid ./TIGR6/ob2sb_TIGR6/sb_by_ob.dist20.synteny_list --orthid ./TIGR6/OrthOrNot/OBRACH_SBICO/OBRACH.ortholog.table --outdir ./TIGR6/SytnteyOrNot --project OBRACH
perl SyntenyOrNot.pl --syntid ./TIGR6/ob2sb_TIGR6/ob_by_sb.dist20.synteny_list --orthid ./TIGR6/OrthOrNot/OBRACH_SBICO/SBICO.ortholog.table --outdir ./TIGR6/SytnteyOrNot --project SBICO

echo "rice sorghum"
perl OrthorNot.pl --typedir ./TIGR6/ClusterGeneType --distance ./TIGR6/Distance/TIGR6_SBICO.distance.txt
perl SyntenyOrNot.pl --syntid ./TIGR6/os2sb_TIGR6/os_by_sb.dist20.synteny_list --orthid ./TIGR6/OrthOrNot/TIGR6_SBICO/SBICO.ortholog.table --outdir ./TIGR6/SytnteyOrNot --project SBICO
perl SyntenyOrNot.pl --syntid ./TIGR6/os2sb_TIGR6/sb_by_os.dist20.synteny_list --orthid ./TIGR6/OrthOrNot/TIGR6_SBICO/TIGR6.ortholog.table --outdir ./TIGR6/SytnteyOrNot --project TIGR6

echo "OBa IRGSP_with_predictgene"
perl OrthorNot.pl --typedir ./IRGSP_withpredictgene/ClusterGeneType --distance ./IRGSP_withpredictgene/Distance/OBRACH_IRGSP.distance.txt > log 2> log2 &
perl SyntenyOrNot.pl --syntid ./IRGSP_withpredictgene/ob2os_IRGSP_withpredictgene/os_by_ob.dist20.synteny_list --orthid ./IRGSP_withpredictgene/OrthOrNot/OBRACH_IRGSP/OBRACH.ortholog.table --outdir ./IRGSP_withpredictgene/ --project OBRACH
perl SyntenyOrNot.pl --syntid ./IRGSP_withpredictgene/ob2os_IRGSP_withpredictgene/ob_by_os.dist20.synteny_list --orthid ./IRGSP_withpredictgene/OrthOrNot/OBRACH_IRGSP/IRGSP.ortholog.table --outdir ./IRGSP_withpredictgene/ --project IRGSP  

echo "get syn gene gff for gramene genebuilder"
perl getidgff.pl --list TIGR6_3way/SyntenyOrNot/OBRACH.syn.ortholog.table --gff /home/jfchen/FFproject/seqlib/BGI_analysis_data/brachyantha/superscaffold/Oryza_brachyantha.genome.super_scaffold.v1.0.glean.gff.150_filter.2K.gff --output OBa.syn.gff

echo "OBa_gramene_v1_final8"
perl ClusterGeneType.pl --cluster ../input/output/all_vs_all.blast.m8.solar.forHC.hcluster --cds ../input/all.cds > log 2> log2 & 
perl ClusterGeneType.pl --cluster ../input/output/all_vs_all.blast.m8.solar.forHC.hcluster --cds ../input/input/all.cds > log 2> log2 &

echo "final6, best4"
perl ClusterGeneType.pl --cluster ../input/result_Gramene_v1_final6_best4/all_vs_all.blast.m8.solar.forHC.hcluster --cds ../input/all.final6.best4.cds > log 2> log2 &

echo "final8_2"
./sumdistance3way.sh > summary
perl ClusterGeneType.pl --cluster ../input/output/all_vs_all.blast.m8.solar.forHC.hcluster --cds ../input/input/all.cds > log 2> log2 &
perl OrthorNot.pl --typedir ./TIGR6_GRAMENE_v1_final8_2/ClusterGeneType --distance TIGR6_GRAMENE_v1_final8_2/Distance/OBRACH_TIGR6.distance.txt > log 2> log2 &
perl SyntenyOrNot.pl --syntid TIGR6_GRAMENE_v1_final8_2/ob2os_TIGR6_Gramene_v1_final8_2/os_by_ob.dist20.synteny_list --orthid TIGR6_GRAMENE_v1_final8_2/OrthOrNot/OBRACH_TIGR6/OBRACH.ortholog.table --outdir TIGR6_GRAMENE_v1_final8_2/ --project OBRACH
perl SyntenyOrNot.pl --syntid TIGR6_GRAMENE_v1_final8_2/ob2os_TIGR6_Gramene_v1_final8_2/ob_by_os.dist20.synteny_list --orthid TIGR6_GRAMENE_v1_final8_2/OrthOrNot/OBRACH_TIGR6/TIGR6.ortholog.table --outdir TIGR6_GRAMENE_v1_final8_2/ --project TIGR6


echo "final8_3"
./sumdistance3way.sh > summary
perl ClusterGeneType.pl --cluster ../input/output/all_vs_all.blast.m8.solar.forHC.hcluster --cds ../input/input/all.cds > log 2> log2 &
perl OrthorNot.pl --typedir ./TIGR6_GRAMENE_v1_final8_3/ClusterGeneType --distance TIGR6_GRAMENE_v1_final8_3/Distance/OBRACH_TIGR6.distance.txt > log 2> log2 &
perl SyntenyOrNot.pl --syntid TIGR6_GRAMENE_v1_final8_3/ob2os_TIGR6_Gramene_v1_final8_3/os_by_ob.dist20.synteny_list --orthid TIGR6_GRAMENE_v1_final8_3/OrthOrNot/OBRACH_TIGR6/OBRACH.ortholog.table --outdir TIGR6_GRAMENE_v1_final8_3/ --project OBRACH
perl SyntenyOrNot.pl --syntid TIGR6_GRAMENE_v1_final8_3/ob2os_TIGR6_Gramene_v1_final8_3/ob_by_os.dist20.synteny_list --orthid TIGR6_GRAMENE_v1_final8_3/OrthOrNot/OBRACH_TIGR6/TIGR6.ortholog.table --outdir TIGR6_GRAMENE_v1_final8_3/ --project TIGR6

echo "final9"
./sumdistance3way.sh > summary
perl ClusterGeneType.pl --cluster ../input/output/all_vs_all.blast.m8.solar.forHC.hcluster --cds ../input/input/all.cds > log 2> log2 &
perl OrthorNot.pl --typedir ./TIGR6_GRAMENE_v1_final9/ClusterGeneType --distance TIGR6_GRAMENE_v1_final9/Distance/OBRACH_TIGR6.distance.txt > log 2> log2 &
perl SyntenyOrNot.pl --syntid TIGR6_GRAMENE_v1_final9/ob2os_TIGR6_Gramene_v1_final9/os_by_ob.dist20.synteny_list --orthid TIGR6_GRAMENE_v1_final9/OrthOrNot/OBRACH_TIGR6/OBRACH.ortholog.table --outdir TIGR6_GRAMENE_v1_final9/ --project OBRACH
perl SyntenyOrNot.pl --syntid TIGR6_GRAMENE_v1_final9/ob2os_TIGR6_Gramene_v1_final9/ob_by_os.dist20.synteny_list --orthid TIGR6_GRAMENE_v1_final9/OrthOrNot/OBRACH_TIGR6/TIGR6.ortholog.table --outdir TIGR6_GRAMENE_v1_final9/ --project TIGR6

echo "final9a"
./sumdistance3way.sh > summary
perl ClusterGeneType.pl --cluster ../input/output/all_vs_all.blast.m8.solar.forHC.hcluster --cds ../input/input/all.cds > log 2> log2 &
perl OrthorNot.pl --typedir ./TIGR6_GRAMENE_v1_final9a/ClusterGeneType --distance TIGR6_GRAMENE_v1_final9a/Distance/OBRACH_TIGR6.distance.txt > log 2> log2 &
perl SyntenyOrNot.pl --syntid TIGR6_GRAMENE_v1_final9a/ob2os_TIGR6_Gramene_v1_final9a/os_by_ob.dist20.synteny_list --orthid TIGR6_GRAMENE_v1_final9a/OrthOrNot/OBRACH_TIGR6/OBRACH.ortholog.table --outdir TIGR6_GRAMENE_v1_final9a/ --project OBRACH
perl SyntenyOrNot.pl --syntid TIGR6_GRAMENE_v1_final9a/ob2os_TIGR6_Gramene_v1_final9a/ob_by_os.dist20.synteny_list --orthid TIGR6_GRAMENE_v1_final9a/OrthOrNot/OBRACH_TIGR6/TIGR6.ortholog.table --outdir TIGR6_GRAMENE_v1_final9a/ --project TIGR6

echo "final9_best"
./sumdistance3way.sh > summary
perl ClusterGeneType.pl --cluster ../input/output/all_vs_all.blast.m8.solar.forHC.hcluster --cds ../input/input/all.cds > log 2> log2 &
perl OrthorNot.pl --typedir ./TIGR6_GRAMENE_v1_final9_best/ClusterGeneType --distance TIGR6_GRAMENE_v1_final9_best/Distance/OBRACH_TIGR6.distance.txt > log 2> log2 &
perl SyntenyOrNot.pl --syntid TIGR6_GRAMENE_v1_final9_best/ob2os_TIGR6_Gramene_v1_final9_best/os_by_ob.dist20.synteny_list --orthid TIGR6_GRAMENE_v1_final9_best/OrthOrNot/OBRACH_TIGR6/OBRACH.ortholog.table --outdir TIGR6_GRAMENE_v1_final9_best/ --project OBRACH
perl SyntenyOrNot.pl --syntid TIGR6_GRAMENE_v1_final9_best/ob2os_TIGR6_Gramene_v1_final9_best/ob_by_os.dist20.synteny_list --orthid TIGR6_GRAMENE_v1_final9_best/OrthOrNot/OBRACH_TIGR6/TIGR6.ortholog.table --outdir TIGR6_GRAMENE_v1_final9_best/ --project TIGR6

echo "get rice gene list that present as syn in final8_2 but not in final8_3"
perl checkid.pl --list TIGR6_GRAMENE_v1_final8_2/SyntenyOrNot/TIGR6.syn.ortholog.table --list2 TIGR6_GRAMENE_v1_final8_3/SyntenyOrNot/TIGR6.syn.ortholog.table > final8_2Overfinal8_3

echo "get related gene that not present in final8_2"
perl getadditionalFFgene.pl --rice final8_2Overfinal8_3 --distance TIGR6_GRAMENE_v1_final8_2/Distance/OBRACH_TIGR6.distance.txt --gff /home/jfchen/FFproject/seqlib/BGI_analysis_data/brachyantha/gramene/input/OBa.gramene.FINAL8_3.bestmodel.gff > OBa.final8_2.additional

echo "get gff of these gene, use as additional model to merge with final8_3_best using besttranscript.pl on 160"
perl getidgff2.pl -l OBa.final8_2.additional -g /home/jfchen/FFproject/seqlib/BGI_analysis_data/brachyantha/gramene/input/OBa.gramene.FINAL8_2.bestmodel.gff -o OBa.final8_2.additional.gff

echo "final8_3_last"
perl ClusterGeneType.pl --cluster ../input/output/all_vs_all.blast.m8.solar.forHC.hcluster --cds ../input/input/all.cds > log 2> log2 &


echo "get shared gff for mcscan"
perl getidgff.pl -l TIGR6_3way/ClusterGeneType/TIGR6.shared.table -g /home/jfchen/FFproject/seqlib/BGI_analysis_data/tigr_rice/tigr.all.final.gff -o TIGR6.shared.gff
perl getidgff.pl -l TIGR6_3way/ClusterGeneType/OBRACH.shared.table -g /home/jfchen/FFproject/FFgenome/06.draw/chralign/dotplot/input/OBa.all.chr.gff -o OBRACH.shared.gff
perl getidgff.pl -l TIGR6_3way/ClusterGeneType/SBICO.shared.table -g /home/jfchen/FFproject/FFgenome/06.draw/chralign/dotplot/input/Sb.final.gff -o SBICO.shared.gff


