echo "test run"
perl PrepareOrthFas.pl --orth colinear.list --seq1 ../input/tigr.all.cds.final.fa --seq2 ../input/Gramene.cds.fa >log 2>log2 &

perl /home/jfchen/FFproject/FFgenome/03.evolution/kaks_pairwise/bin/script_KaKs_calculator/pipe_kaks.pl ./colinear > KaKs.summary

echo "collinear gene run"
perl PrepareOrthFas.pl --orth ../input/chr01.list --seq1 ../input/tigr.all.cds.final.fa --seq2 ../input/Gramene.cds.fa >log 2>log2 &
perl /home/jfchen/FFproject/FFgenome/03.evolution/kaks_pairwise/bin/script_KaKs_calculator/pipe_kaks.pl ./chr01 > chr01.kaks.summary 2> log2 &

echo "non collinear 2 Os pep best hit"
perl PrepareOrthFas.pl -orth rice.Other.noncollinear.nc2Os.solar.best.identity -seq1 /data/jfchen/Final_work/seqlib/tigr6.1/Os.cds -seq2 /data/jfchen/Final_work/seqlib/tigr6.1/Os.cds > log 2> log2 &

perl /home/jfchen/FFproject/FFgenome/03.evolution/kaks_pairwise/bin/script_KaKs_calculator/pipe_kaks.pl ./rice.Other.noncollinear.nc2Os.solar.best > rice.Other.noncollinear.nc2Os.solar.best.kaks.summary 2> log2 &




