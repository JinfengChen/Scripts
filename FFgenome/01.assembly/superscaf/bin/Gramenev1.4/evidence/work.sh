perl /home/jfchen/FFproject/tools/bin/runblastall_multicpu.pl -p blastp -i ../Gramene.pep.fa -d plant_protein.fa > log 2> log2 &

perl /home/jfchen/FFproject/tools/bin/runblastall_multicpu.pl -p blastn -i ../Gramene.cds.fa -d cuflink.tran.fa > log 2> log2 &

perl evidence.pl --gene ../Gramene.cds.fa --protein Gramene.pep.fa.blasttable --transcript Gramene.cds.fa.blasttable --iprscan Gramene.final.v1.4.pep.iprscan > log 2> log2 &



