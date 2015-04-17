perl runblastall_multicpu_m8.pl -p blastn -i ../input/brachyantha.bes -d /home/jfchen/FFproject/seqlib/BGI_analysis_data/tigr_rice/all.con > log 2 >log2 &

perl statm8.pl --inf ../input/infor_FPC.txt --m8 brachyantha.bes.blast.m8

                                                                                                       
blat ../input/Gramene.chr.fa ../input/brachyantha.bes brachyantha.bes2FF.psl -noHead -fastMap > log 2> log2 &
perl insertSD.pl --psl brachyantha.bes2FF.psl --inf ../input/infor_FPC.txt > log 2> log2 &

perl runblastall_multicpu.pl -p blastn -i ../input/brachyantha.bes -d /home/jfchen/FFproject/seqlib/BGI_analysis_data/tigr_rice/all.con > log 2> log2 &

perl classifyclone.pl --blasttable test.blasttable --inf ../input/infor_FPC.txt

perl filterblasttable.pl --blasttable brachyantha.bes.blasttable --inf ../input/infor_FPC.txt > filter.log &
perl classifyclone.pl --inf ../input/infor_FPC.txt --blasttable brachyantha.bes.filter.blasttable > clone.log 2> clone.log2 &


perl m8tblasttable.pl --query ../input/brachyantha.bes --target ../input/Gramene.chr.fa --blastm8 test.blast.m8

perl clusterbespair.pl --classify brachyantha.bes.filter.classify > log 2> log2 &

blat ../input/all.con ../input/glaberrima.bes glaberrima.bes.blast.m8 -out=blast8 -minScore=160 -noHead -fastMap > log 2> log2 &

echo "blast BES to rice genome"
blat ../input/all.con ../input/glaberrima.bes glaberrima.bes.blast.m8 -out=blast8 -minScore=160 -noHead -fastMap > log 2> log2 &
perl runblastall_multicpu_m8.pl -p blastn -i ../input/punctata.bes -d /home/jfchen/FFproject/seqlib/BGI_analysis_data/tigr_rice/all.con > blast.log 2> blast.log2 &

echo "m8 to blasttable and filter"
perl m8tblasttable.pl --query ../input/glaberrima.bes --target ../input/all.con --blastm8 glaberrima.bes.blast.m8
perl filterblasttable.pl --blasttable glaberrima.bes.blasttable --identity 0.95 --coverage 0.5


echo "classify and cluster"
perl classifyclone.pl --inf ../input/infor_FPC.txt --blasttable glaberrima.bes.filter.blasttable > log 2> log2 &
perl clusterbespairV2.pl --classify glaberrima.bes.filter.classify > log 2> log2 &

echo "analyze cluster"
perl analyzacluster.pl --cluster glaberrima.bes.filter.cluster --gff ../input/tigr.all.final.gff

