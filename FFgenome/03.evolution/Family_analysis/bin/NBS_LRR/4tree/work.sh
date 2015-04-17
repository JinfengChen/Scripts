cp ../OBa/*.fa ./
cp ../tigr6/*.fa ./
cp ../SB/*.fa ./
./merge.sh
rm SB.* OBa.* rice.*
perl ../../DrawTree.pl --protein CC-NBS.fa --align > log 2> log2 &
perl ../../DrawTree.pl --protein NBS.fa --align > log3 2> log4 &
perl ../../DrawTree.pl --protein CC-NBS-LRR.fa --align > log 2> log2 &
perl ../../DrawTree.pl --protein NBS-LRR.fa --align > log 2> log2 &

perl ../../DrawTree.pl --protein NBS.rice.fa --align --nj > log 2> log2 &

hmmalign --trim --outformat PSIBLAST ../../../input/hmm2/NB-ARC.hmm NBS.rice.fa > log 2> log2 &
perl ../../hmmalign2fasta.pl --hmmalign NBS.rice.fa.hmmalign > NBS.rice.fa.hmmalign.fasta 2> log2 &
/home/jfchen/FFproject/tools/FastTree/FastTree NBS.rice.fa.hmmalign.fasta > NBS.rice.fa.hmmalign.fasta.tree 2> log2 &



