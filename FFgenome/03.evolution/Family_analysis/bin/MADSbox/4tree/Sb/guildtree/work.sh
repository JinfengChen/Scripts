perl /home/jfchen/software/tree/tree.pl --leaf SB.MADS.M1.tree.nwk | grep Sb > SB.MADS.M1.ID
perl /home/jfchen/software/tree/tree.pl --leaf SB.MADS.M2.tree.nwk | grep Sb > SB.MADS.M2.ID
perl /home/jfchen/software/tree/tree.pl --leaf SB.MADS.M3.tree.nwk | grep Sb > SB.MADS.M3.ID
perl /home/jfchen/software/tree/tree.pl --leaf SB.MADS.MIKC.tree.nwk | grep Sb > SB.MADS.MIKC.ID

perl ../../../../getidseq.pl --list SB.MADS.M1.ID --fasta ../SB.MADS.4tree.fa --output SB.MADS.M1.fa
perl ../../../../getidseq.pl --list SB.MADS.M2.ID --fasta ../SB.MADS.4tree.fa --output SB.MADS.M2.fa
perl ../../../../getidseq.pl --list SB.MADS.M3.ID --fasta ../SB.MADS.4tree.fa --output SB.MADS.M3.fa
perl ../../../../getidseq.pl --list SB.MADS.MIKC.ID --fasta ../SB.MADS.4tree.fa --output SB.MADS.MIKC.fa

