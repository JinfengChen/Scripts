perl ../../DrawTree.pl --protein Other.fa --align > log 2> log2 &
perl ../../DrawTree.pl --protein Kelch.fa --align > log 2> log2 &
perl ../../DrawTree.pl --protein FBA.fa --align > log 2> log2 &
perl ../../DrawTree.pl --protein DUF295.fa --align > log 2> log2 &
perl ../../DrawTree.pl --protein LRR.fa --align > log 2> log2 &
perl ../../DrawTree.pl --protein Unknown.fa --align > log 2> log2 &
hmmalign --trim --outformat PSIBLAST ../../../input/hmm/F-box.hmm Unknown.fa > Unknown.fa.hmmalign
perl hmmalign2phy.pl --hmmalign Unknown.fa.hmmalign > Unknown.fa.hmmalign.phy
perl ../../DrawTree.pl --protein Unknown.fa.hmmalign.phy

