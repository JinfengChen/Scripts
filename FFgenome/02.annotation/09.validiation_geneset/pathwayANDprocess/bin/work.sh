perl findOrthGene.pl --geneid ../input/RiceGene.id --ortholog ../input/OBRACH_TIGR6.distance.txt --solor ../input/all_vs_all.blast.m8.solar --hcluster ../input/all_vs_all.blast.m8.solar.forHC.hcluster --gff ../input/tigr.all.final.gff --project RiceGene > RiceGene.miss 2> log2 &

perl findOrthGene.pl --geneid ../input/StarchPathway.id --ortholog ../input/OBRACH_TIGR6.distance.txt --solor ../input/all_vs_all.blast.m8.solar --hcluster ../input/all_vs_all.blast.m8.solar.forHC.hcluster --gff ../input/tigr.all.final.gff --project Starch > Starch.miss 2> log2 &

perl findOrthGene.pl --geneid ../input/flowering.id --ortholog ../input/OBRACH_TIGR6.distance.txt --solor ../input/all_vs_all.blast.m8.solar --hcluster ../input/all_vs_all.blast.m8.solar.forHC.hcluster --gff ../input/tigr.all.final.gff --project Flowering > Flowering.miss 2> log2 &


perl findOrthGene.pl --geneid ../input/RiceGene.id --ortholog ../input/OBRACH_TIGR6.distance.txt --solor ../input/result_TIGR6_3way/all_vs_all.blast.m8.solar --hcluster ../input/result_TIGR6_3way/all_vs_all.blast.m8.solar.forHC.hcluster --gff ../input/tigr.all.final.gff --project RiceGene > RiceGene.miss
perl findOrthGene.pl --geneid ../input/flowering.id --ortholog ../input/OBRACH_TIGR6.distance.txt --solor ../input/result_TIGR6_3way/all_vs_all.blast.m8.solar --hcluster ../input/result_TIGR6_3way/all_vs_all.blast.m8.solar.forHC.hcluster --gff ../input/tigr.all.final.gff --project Flowering > Flowering.miss
perl findOrthGene.pl --geneid ../input/StarchPathway.id --ortholog ../input/OBRACH_TIGR6.distance.txt --solor ../input/result_TIGR6_3way/all_vs_all.blast.m8.solar --hcluster ../input/result_TIGR6_3way/all_vs_all.blast.m8.solar.forHC.hcluster --gff ../input/tigr.all.final.gff --project Starch > Starch.miss
