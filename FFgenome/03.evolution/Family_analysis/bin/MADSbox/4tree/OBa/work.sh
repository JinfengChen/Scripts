cat ../../OB.MADS.fa ../../OS.MADS.representive.fa > OBa.MADS.4tree.fa
perl ../../../DrawTree.pl --protein OBa.MADS.4tree.fa --align > log 2> log2 &

