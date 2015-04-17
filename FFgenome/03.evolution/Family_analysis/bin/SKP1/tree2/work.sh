perl /home/jfchen/software/tree/tree.pl --leaf tree2.tre | grep "LOC" > OS.tree2.leaf
perl /home/jfchen/software/tree/tree.pl --leaf tree2.tre | grep "OB" > OB.tree2.leaf
perl ../../getidseq.pl --list OS.tree2.leaf --fasta ../skp1.fa --out OS.tree2.fa
perl ../../getidseq.pl --list OB.tree2.leaf --fasta ../skp1.fa --out OB.tree2.fa
cat OS.tree2.fa OB.tree2.fa > tree2.fa
perl ../../DrawTree.pl --protein tree2.fa --align > log 2> log2 &

