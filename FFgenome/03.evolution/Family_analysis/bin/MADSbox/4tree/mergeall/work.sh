cp ../OBa/guidetree/*.fa ./
cp ../Sb/guildtree/*.fa ./
cp ../rice/type/*.fa ./
cat *.M1.* > MADS.M1.fa
cat *.M2.* > MADS.M2.fa
cat *.M3.* > MADS.M3.fa
cat *.MIKC.* > MADS.MIKC.fa

perl ../../../DrawTree.pl --protein MADS.M1.fa --align > log 2> log2 &

perl ../../../DrawTree.pl --protein MADS.M2.fa --align > log 2> log2 &
perl ../../../DrawTree.pl --protein MADS.M3.fa --align > log 2> log2 &
