cat ../OBa.Blast.found.fa ../../rice_RLK_LRR/RLK_LRR_rice.representive.fa > OBa.RLK-LRR.4tree.fa
hmmalign --trim --outformat PSIBLAST ../../../../input/Pfam/Pkinase.hmm OBa.RLK-LRR.4tree.fa > OBa.RLK-LRR.4tree.fa.hmmalign
perl ../../../hmmalign2phy.pl --hmmalign OBa.RLK-LRR.4tree.fa.hmmalign > OBa.RLK-LRR.4tree.fa.hmmalign.phy
perl ../../../DrawTree.pl --protein OBa.RLK-LRR.4tree.fa.hmmalign.phy > log 2> log2 &

perl ../../../hmmalign2fasta.pl --hmmalign OBa.RLK-LRR.4tree.fa.hmmalign > OBa.RLK-LRR.4tree.fa.hmmalign.fa
treebest nj OBa.RLK-LRR.4tree.fa.hmmalign.fa > nj.tree
