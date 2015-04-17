perl ../bin/convertBLAST.pl --blast 4way.blasttable --project 4way &

echo "convert GFF"
perl ../bin/convertGFF.pl --gff Sb.final.gff --project Sb
perl ../bin/convertGFF.pl --gff OBa.all.chr.gff --project Ob
perl ../bin/convertGFF.pl --gff tigr.all.final.gff --project Os
echo "convert blastm8 from hcluster pipeline"
perl ../bin/convertBLAST4m8.pl --blastm8 all_vs_all.blast.m8 --project 3way

echo "run"
cp 3way.* ../bin/
cd ../bin
perl runMcscan.pl -mcl --project 3way > log 2> log2 &


echo "2way"
grep -v "Sb" 3way.blast > 2way.blast
grep -v "Sb" 3way.gff > 2way.gff

echo "gramene_v1_final2"
perl ../bin/convertBLAST4m8.pl --blastm8 gramene_v1_final2.blast.m8 --project final2
perl ../bin/convertGFF.pl --gff Gramene.chr.gff --project Gr
cat Gr.gff Os.gff > final2.gff
grep -v "Sb" final2.blast > temp.blast
mv temp.blast final2.blast


echo "gramene_v1_final6"
perl ../bin/convertBLAST4m8.pl --blastm8 gramene_v1.final6.blastm8 --project final6
perl ../bin/convertGFF.pl --gff Gramene.chr.gff --project Gr
cat Gr.gff Os.gff > final6.gff
grep -v "Sb" final6.blast > temp.blast
mv temp.blast final6.blast

cp final6.* ../bin/
cd ../bin
perl runMcscan.pl --mcl --project final6 > log 2> log2 &
perl sumBlock.pl --block final6.blocks > final6.blocks.summary

echo "gramene_v1_final6_best2"
perl ../bin/convertBLAST4m8.pl --blastm8 gramene.best2.blastm8 --project best2
perl ../bin/convertGFF.pl --gff Gramene.chr.gff --project Gr
cat Gr.gff Os.gff > best2.gff

grep -v "Sb" best2.blast > temp.blast
mv temp.blast best2.blast
cp best2.* ../bin/
cd ../bin/
perl runMcscan.pl --mcl --project best2 > log 2> log2 &
perl sumBlock.pl --block best2.blocks > best2.blocks.summary

perl check.pl --list2 gramene_v1_final6_best1/best1.blocks.table --list gramene_v1_final6_best2/best2.blocks.table > best2overbest1.sum 
perl check.pl --list gramene_v1_final6_best1/best1.blocks.table --list2 gramene_v1_final6_best2/best2.blocks.table > best1overbest2.sum
perl check.pl --list 2way_chr/2way.blocks.table --list2 gramene_v1_final6_best2/best2.blocks.table > BGIoverbest2.sum 
perl check.pl --list2 2way_chr/2way.blocks.table --list gramene_v1_final6_best2/best2.blocks.table > best2overBGI.sum

echo "gramene_v1_final6_best3"
perl ../bin/convertBLAST4m8.pl --blastm8 gramene.best3.blastm8 --project best3
perl ../bin/convertGFF.pl --gff Gramene.chr.gff --project Gr
cat Gr.gff Os.gff > best3.gff

grep -v "Sb" best3.blast > temp.blast
mv temp.blast best3.blast
cp best3.* ../bin/
cd ../bin/
perl runMcscan.pl --mcl --project best3 > log 2> log2 &
perl sumBlock.pl --block best3.blocks > best3.blocks.summary

echo "format all_vs_all.blastm8 to normal name format and exclude species unlike"
perl ../bin/BLASTm8name.pl --blast all_vs_all.blast.m8 --exclude SBICO > log 2> log2 &

echo "prepare bed file for chromosome"
perl ../bin/convertGFF.pl --gff tigr.all.final.gff --project Os
perl ../bin/convertGFF.pl --gff OBa.all.chr.gff --project Ob
perl ../bin/convertGFF.pl --gff Sb.final.gff --project Sb
perl ../bin/convertGFF.pl --gff Gramene.FINAL8_3_last.chr.v1.gff --project Gr
perl ../bin/convertGFF.pl --gff Gramene.chr.gff --project Ob
perl ../bin/splitGFF.pl Os.bed 
perl ../bin/splitGFF.pl Ob.bed 
perl ../bin/splitGFF.pl Sb.bed
perl ../bin/splitGFF.pl Gr.bed

echo "gramene_final8_3_last, TIGR6, sorghum"
perl ../bin/convertBLAST4m8.pl --blastm8 all_vs_all.blast.m8 --project 3way > log 2> log2 &


echo "gramenev1.4, tigr6, sorghum"
perl ../bin/convertBLAST4m8.pl --blastm8 all_vs_all.blast.m8 --project 3wayv1.4 >log 2>log2

