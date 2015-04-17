perl ../../../DrawTree.pl --protein Other.fa --align > log 2> log2 &
perl ../../../DrawTree.pl --protein DUF295.fa --align > log 2> log2 &
perl ../../../DrawTree.pl --protein LRR-FDB.fa --align > log 2> log2 &
perl ../../../DrawTree.pl --protein FBA-Kelch.fa --align > log 2> log2 &

hmmalign --trim --outformat PSIBLAST ../../../../input/hmm/F-box.hmm Unknown.fa > Unknown.fa.hmmalign
perl ../../../hmmalign2phy.pl --hmmalign Unknown.fa.hmmalign > Unknown.fa.hmmalign.phy


treebest nj DUF295.fa.muscle > nj.log 2> nj.log2 &
perl /home/jfchen/software/tree/tree_plot.pl --png nj.log

