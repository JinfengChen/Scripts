echo "test example"
perl synteny_table.pl ../input/os_sb.ort > log 2> log2 &

echo "run ob_os" 
perl synteny_table.pl ../input/os_ob.ort > log 2> log2 &

echo "OBa and TIGR6"
perl inf2ort.pl --inf /home/jfchen/FFproject/FFgenome/03.evolution/genefamily_hcluster/sumType/bin/TIGR6/Distance/OBRACH_TIGR6.distance.txt --gff1 ../input/OBa.all.chr.gff --gff2 ../input/tigr.all.final.gff --project ob_os
perl synteny_table.pl ../input/ob_os.ort > log 2> log2 &

echo "TIGR6 and SB"
perl inf2ort.pl --inf /home/jfchen/FFproject/FFgenome/03.evolution/genefamily_hcluster/sumType/bin/TIGR6/Distance/TIGR6_SBICO.distance.txt --gff1 ../input/tigr.all.final.gff --gff2 ../input/Sb.final.gff --project os_sb
perl synteny_table.pl ../input/os_sb.ort > log 2> log2 &

echo "OBa and SB"
perl inf2ort.pl --inf /home/jfchen/FFproject/FFgenome/03.evolution/genefamily_hcluster/sumType/bin/TIGR6/Distance/OBRACH_SBICO.distance.txt --gff1 ../input/OBa.all.chr.gff --gff2 ../input/Sb.final.gff --project ob_sb
perl synteny_table.pl ../input/ob_sb.ort > log 2> log2 &

echo "OBa and IRGSPwithpredictgene"
perl inf2ort.pl --inf /home/jfchen/FFproject/FFgenome/03.evolution/genefamily_hcluster/sumType/bin/IRGSP_withpredictgene/Distance/OBRACH_IRGSP.distance.txt --gff1 ../input/OBa.all.chr.gff --gff2 ../input/IRGSP.all.gff --project ob_os
perl synteny_table.pl ../input/ob_os.ort > log 2> log2 &

echo "OBa and TIGR6 for drawACT"
perl synteny_table.pl ../input/ob_os.ort > log 2> log2 &

echo "OBa and TIGR6, final8_3"
perl inf2ort.pl --inf /home/jfchen/FFproject/FFgenome/03.evolution/genefamily_hcluster/sumType/bin/TIGR6_GRAMENE_v1_final8_3/Distance/OBRACH_TIGR6.distance.txt --gff1 ../input/Gramene.FINAL8_3.chr.gff --gff2 ../input/tigr.all.final.gff --project ob_os > log 2> log2 &
perl synteny_table.pl ../input/ob_os.ort > log 2> log2 &

echo "OBa and TIGR6, final9"
perl inf2ort.pl --inf /home/jfchen/FFproject/FFgenome/03.evolution/genefamily_hcluster/sumType/bin/TIGR6_GRAMENE_v1_final9/Distance/OBRACH_TIGR6.distance.txt --gff1 ../input/Gramene.FINAL9.chr.gff --gff2 ../input/tigr.all.final.gff --project ob_os > log 2> log2 &
perl synteny_table.pl ../input/ob_os.ort > log 2> log2 &

echo "OBa and TIGR6, final8_2"
perl inf2ort.pl --inf /home/jfchen/FFproject/FFgenome/03.evolution/genefamily_hcluster/sumType/bin/TIGR6_GRAMENE_v1_final8_2/Distance/OBRACH_TIGR6.distance.txt --gff1 ../input/Gramene.FINAL8_2.chr.gff --gff2 ../input/tigr.all.final.gff --project ob_os > log 2> log2 &
perl synteny_table.pl ../input/ob_os.ort > log 2> log2 &

echo "OBa and TIGR6, final9a"
perl inf2ort.pl --inf /home/jfchen/FFproject/FFgenome/03.evolution/genefamily_hcluster/sumType/bin/TIGR6_GRAMENE_v1_final9a/Distance/OBRACH_TIGR6.distance.txt --gff1 ../input/Gramene.FINAL9a.chr.gff --gff2 ../input/tigr.all.final.gff --project ob_os > log 2> log2 &
perl synteny_table.pl ../input/ob_os.ort > log 2> log2 &

echo "OBa and TIGR6, final9_best"
perl inf2ort.pl --inf /home/jfchen/FFproject/FFgenome/03.evolution/genefamily_hcluster/sumType/bin/TIGR6_GRAMENE_v1_final9_best/Distance/OBRACH_TIGR6.distance.txt --gff1 ../input/Gramene.FINAL9_best.chr.gff --gff2 ../input/tigr.all.final.gff --project ob_os > log 2> log2 &
perl synteny_table.pl ../input/ob_os.ort > log 2> log2 &

echo "OBa and TIGR6, final8_3_last"
perl inf2ort.pl --inf /home/jfchen/FFproject/FFgenome/03.evolution/genefamily_hcluster/sumType/bin/TIGR6_GRAMENE_v1_final8_3_last/Distance/OBRACH_TIGR6.distance.txt --gff1 ../input/Gramene.FINAL8_3_last.chr.gff --gff2 ../input/tigr.all.final.gff --project ob_os > log2 > log2 &
perl synteny_table.pl ../input/ob_os.ort > log 2> log2 &

