cat ../SB.Blast.found.fa ../../rice_RLK_LRR/RLK_LRR_rice.representive.fa > SB.RLK-LRR.4tree.fa
hmmalign --trim --outformat PSIBLAST ../../../../input/Pfam/Pkinase.hmm SB.RLK-LRR.4tree.fa > SB.RLK-LRR.4tree.fa.hmmalign
perl ../../../hmmalign2phy.pl --hmmalign SB.RLK-LRR.4tree.fa.hmmalign > SB.RLK-LRR.4tree.fa.hmmalign.phy
perl ../../../DrawTree.pl --protein SB.RLK-LRR.4tree.fa.hmmalign.phy > log 2> log2 &

