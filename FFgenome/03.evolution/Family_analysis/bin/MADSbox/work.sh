echo "MADSbox contain SRF-TF/F00319 pfam domain"

perl ../Get2PfamGeneID4hmmer3.pl --hmmer3 ../tigr.pfam3.hmmer3 --domain1 PF00319 --domain2 PF00069 > MADS.rice
perl ../Get2PfamGeneID4hmmer3.pl --hmmer3 ../OBa.pfam3.hmmer3 --domain1 PF00319 --domain2 PF00069 > MADS.OBa
perl ../Get2PfamGeneID4hmmer3.pl --hmmer3 ../Sb.pfam3.hmmer3 --domain1 PF00319 --domain2 PF00069 > MADS.SB
perl ../getidseq.pl --list MADS.rice --fasta ../../input/tigr.all.pep.final.fa -o OS.MADS.fa
perl ../getidseq.pl --list MADS.OBa --fasta ../../input/Oryza_brachyantha.genome.super_scaffold.v1.0.glean.gff.150_filter.2K.pep -o OB.MADS.fa
perl ../getidseq.pl --list MADS.SB --fasta ../../input/Sb.final.pep -o SB.MADS.fa

perl ../DrawTree.pl --protein MADS.fa --align > log 2> log2 &

perl checkID.pl --id4pfam MADS.rice --id4pub ./MADS_BMC2007/TIGR.MADS.BMC2007.ID
perl AssignType.pl --fasta OS.MADS.fa --id MADS_BMC2007/TIGR.MADS.BMC2007.txt --output OS.MADS.ID.fa
perl ../getidseq.pl --list OS.MADS.representive.id --fasta OS.MADS.fa --o OS.MADS.representive.fa

