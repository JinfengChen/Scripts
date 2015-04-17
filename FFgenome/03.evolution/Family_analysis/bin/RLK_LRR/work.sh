echo "recept like kanise contain Pkinase/PF00069 and LRR/PF00560 domain"
perl ../Get2PfamGeneID4hmmer3.pl --hmmer3 ../tigr.pfam3.hmmer3 --domain1 PF00069 --domain2 PF00560 > Pkinase.OS.ID
perl ../getidseq.pl --list Pkinase.OS.ID --fasta ../../input/tigr.all.pep.final.fa --output Pkinase.OS.fa

echo "we use blast to retrive homolous sequence with similarity > 60%, by checking in rice we miss only one and found 8 additional copy, which could be kick out when building tree"
perl FindRLKLRR.pl --fasta Pkinase.OS.fa --refer ./rice_RLK_LRR/Rice.RLK-LRR.fa --project rice
perl checkID.pl --id4pfam Blast.found -id4pub ./rice_RLK_LRR/Rice.RLK-LRR.ID

perl /home/jfchen/FFproject/tools/bin/fastaDeal.pl --attr id rice.RLK-LRR.fa > rice.RLK-LRR.ID
perl checkID.pl --id4pfam rice.RLK-LRR.ID -id4pub ./rice_RLK_LRR/Rice.RLK-LRR.ID 
perl checkID.pl --id4pfam PF00560.PF00069.tigr.pfam3.id -id4pub ./rice_RLK_LRR/Rice.RLK-LRR.ID

echo "OBa"

perl ../Get2PfamGeneID4hmmer3.pl --hmmer3 ../OBa.pfam3.hmmer3 --domain1 PF00069 --domain2 PF00560 > Pkinase.OB.ID
perl ../getidseq.pl --list Pkinase.OB.ID --fasta ../../input/Oryza_brachyantha.genome.super_scaffold.v1.0.glean.gff.150_filter.2K.pep --output Pkinase.OB.fa
perl FindRLKLRR.pl --fasta Pkinase.OB.fa --refer ./rice_RLK_LRR/Rice.RLK-LRR.fa --project OBa
perl FindRLKLRR.pl --fasta Pkinase.OB.fa --refer ./rice_RLK_LRR/Rice.RLK-LRR.fa --type ./rice_RLK_LRR/Rice.RLK-LRR.type --project OBa


echo "glaberrima"
hmmsearch --tblout OGA.pfam3.hmmer3 ../input/Pfam/Pfam-A.hmm ../input/glaberrima.final.pep.fa > log 2> log2 &
perl ../Get2PfamGeneID4hmmer3.pl --hmmer3 ../OGA.pfam3.hmmer3 --domain1 PF00069 --domain2 PF00560 > Pkinase.OGA.ID
perl ../getidseq.pl --list Pkinase.OGA.ID --fasta ../../input/glaberrima.final.pep.fa --output Pkinase.OGA.fa
perl FindRLKLRR.pl --fasta Pkinase.OGA.fa --refer ./rice_RLK_LRR/Rice.RLK-LRR.fa --type ./rice_RLK_LRR/Rice.RLK-LRR.type --project OGA > log 2> log2 &


