perl /home/jfchen/software/tree/tree.pl --leaf OBa.MADS.M1.tree.nwk | grep OB > OBa.MADS.M1.ID
perl /home/jfchen/software/tree/tree.pl --leaf OBa.MADS.M2.tree.nwk | grep OB > OBa.MADS.M2.ID
perl /home/jfchen/software/tree/tree.pl --leaf OBa.MADS.M3.tree.nwk | grep OB > OBa.MADS.M3.ID
perl /home/jfchen/software/tree/tree.pl --leaf OBa.MADS.MIKC.tree.nwk | grep OB > OBa.MADS.MIKC.ID

perl ../../../../getidseq.pl --list OBa.MADS.M1.ID --fasta ../OBa.MADS.4tree.fa --output OBa.MADS.M1.fa
perl ../../../../getidseq.pl --list OBa.MADS.M2.ID --fasta ../OBa.MADS.4tree.fa --output OBa.MADS.M2.fa
perl ../../../../getidseq.pl --list OBa.MADS.M3.ID --fasta ../OBa.MADS.4tree.fa --output OBa.MADS.M3.fa
perl ../../../../getidseq.pl --list OBa.MADS.MIKC.ID --fasta ../OBa.MADS.4tree.fa --output OBa.MADS.MIKC.fa
