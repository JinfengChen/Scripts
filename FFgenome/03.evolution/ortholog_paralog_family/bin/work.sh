echo "test"
perl ortholog_paralog_family_159.pl ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species HUMAN --distance > log 2> log2 &

echo "run"
cp /home/jfchen/FFproject/seqlib/BGI_analysis_data/cds4ortholog_paralog_pipeline/bin/all.cds ./
perl cds2aa.pl ../input/all.cds > all.pep
perl ortholog_paralog_family_159.pl ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species OB --distance > log 2> log2 &

echo "add Poptr and run"
perl cds2aa.pl ../input/all.cds > all.pep
echo "we have all_vs_all.blast.m8.solar.forHC.hcluster from BGI, so skip step 1 and 2."
perl ortholog_paralog_family_159.pl ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species OBRACH --distance --step 34567 > log 2> log2 &

echo "use tigr6 data replace IRGSP"
perl cds2aa.pl ../input/all.cds > all.pep
perl ortholog_paralog_family_159.pl ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species OBRACH --distance > log 2> log2 &

echo "add predict gene into IRGSP gene set, totally 43230 genes"
perl cds2aa.pl ../input/all.cds > ../input/all.pep
perl ortholog_paralog_family_159.pl ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species OBRACH --distance > log 2> log2 &

echo "use only sorghum Tigr6.1 and OBa to find conlinear gene between rice and FF"
perl cds2aa.pl ../input/all.cds > ../input/all.pep
perl ortholog_paralog_family_159.pl --align_rate 0.2 ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species OBRACH --distance > log 2> log2 &

echo "gramene OBa and TIGR6,Sorghum"
perl cds2aa.pl ../input/all.cds > ../input/all.pep
perl ortholog_paralog_family_159.pl --align_rate 0.2 ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species OBRACH --distance --step 1234567 > log 2> log2 &

echo "gramene OBa_gramene_v1_final1 and TIGR6,Sorghum"
perl cds2aa.pl ../input/all.cds > ../input/all.pep
perl ortholog_paralog_family_159.pl --align_rate 0.2 ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species OBRACH --distance --step 1234567 > log 2> log2 &


echo "gramene OBa_gramene_v1_final2 and TIGR6,Sorghum"
perl cds2aa.pl ../input/all.cds > ../input/all.pep
perl ortholog_paralog_family_159.pl --align_rate 0.2 ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species OBRACH --distance --step 1234567 > log 2> log2 &

echo "gramene OBa_gramene_v1_final6 and TIGR6,Sorghum"
perl cds2aa.pl ../input/all.cds > ../input/all.pep
perl ortholog_paralog_family_159.pl --align_rate 0.2 ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species OBRACH --distance --step 1234567 > log 2> log2 &

echo "gramene OBa_gramene-v1_final6 longest transcript fixed but"
perl cds2aa.pl ../input/all.cds > ../input/all.pep
perl ortholog_paralog_family_159.pl --align_rate 0.2 ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species OBRACH --distance --step 2 > log 2> log2 &

echo "gramene OBa_gramene-v1_final6 all transcript fixed but"
perl cds2aa.pl ../input/all.cds > ../input/all.pep
perl ortholog_paralog_family_159.pl --align_rate 0.2 ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species OBRACH --distance --step 2 > log 2> log2 &

echo "gramene v1 final6, best3"
perl ortholog_paralog_family_159.pl --align_rate 0.2 ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species OBRACH --distance --step 2 > log 2> log2 &
perl ortholog_paralog_family_159.pl --align_rate 0.2 ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species OBRACH --distance --step 34567 > log 2> log2 &

echo "gramene v1 final6 best4, adding fgenesh"
perl ortholog_paralog_family_159.pl --align_rate 0.2 ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species OBRACH --distance --step 234567 > log 2> log2 &

echo "gramene v1 final8, adding fgenesh"
perl ortholog_paralog_family_159.pl --align_rate 0.2 ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species OBRACH --distance --step 234567 > log 2> log2 &

echo "gramene v1 final8_2, adding fgenesh and glean"
perl ortholog_paralog_family_159.pl --align_rate 0.2 ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species OBRACH --distance --step 234567 > log 2> log2 &

echo "gramene v1 final8_3,adding fgenesh and glean, non-overlap"
perl ortholog_paralog_family_159.pl --align_rate 0.2 ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species OBRACH --distance --step 234567 > log 2> log2 &

echo "gramene v1 final9, longest not best,adding fgenesh and glean, non-overlap"
perl ortholog_paralog_family_159.pl --align_rate 0.2 ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species OBRACH --distance --step 234567 > log 2> log2 &

echo "gramene v1 final9a, adding fgenesh and glean, non-overlap"
perl ortholog_paralog_family_159.pl --align_rate 0.2 ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species OBRACH --distance --step 23456 > log 2> log2 &

echo "gramene v1 final9_best, adding fgenesh and glean, non-overlap"
perl ortholog_paralog_family_159.pl --align_rate 0.2 ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species OBRACH --distance --step 23456 > log 2> log2 &

echo "gramene v1 final8_3_last,adding fgenesh and glean, non-overlap"
perl ortholog_paralog_family_159.pl --align_rate 0.2 ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species OBRACH --distance --step 23456 > log 2> log2 &

echo "gramenev1.4"
perl ortholog_paralog_family_159.pl --align_rate 0.2 ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species OBRACH --distance --step 1 > log 2> log2 &

perl ortholog_paralog_family_159.pl --align_rate 0.2 ../input/all.pep ../input/all.cds ../input/all.structure ../input/all.annotation ../input/category.txt ../input/species.nh --key_species OJ --distance --step 2 > log 2> log2 &


