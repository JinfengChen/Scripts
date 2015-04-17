echo "SKP1 gene have two conserved domain Skp1_POZ/PF03931 and Skp1/PF01466"
perl ../Get2PfamGeneID.pl --iprscan ../../input/iprscan/OS.iprscan --domain1 PF01466 --domain2 PF03931
perl ../Get2PfamGeneID4hmmer3.pl --hmmer3 ../tigr.pfam3.hmmer3 --domain1 PF01466 --domain2 PF03931

echo "kong skp1 sequence and tree"
perl ../getidseq.pl --list kong.skp1.rice --fasta ../../input/tigr.all.pep.final.fa --o kong.skp1.fa
perl ../Get2PfamGeneID4hmmer3.pl --hmmer3 ../OBa.pfam3.hmmer3 --domain1 PF01466 --domain2 PF03931
perl ../Get2PfamGeneID4hmmer3.pl --hmmer3 ../Sb.pfam3.hmmer3 --domain1 PF01466 --domain2 PF03931
perl ../Get2PfamGeneID4hmmer3.pl --hmmer3 ../tigr.pfam3.hmmer3 --domain1 PF01466 --domain2 PF03931

perl ../getidseq.pl --list PF01466.PF03931.OBa.pfam3.id --fasta ../../input/Oryza_brachyantha.genome.super_scaffold.v1.0.glean.gff.150_filter.2K.pep -o OBa.skp1.fa
perl ../getidseq.pl --list PF01466.PF03931.Sb.pfam3.id --fasta ../../input/Sb.final.pep -o SB.skp1.fa
perl ../getidseq.pl --list PF01466.PF03931.tigr.pfam3.id --fasta ../../input/tigr.all.pep.final.fa -o OS.skp1.fa
perl ../getidseq.pl --list kong.skp1.rice --fasta ../../input/tigr.all.pep.final.fa --o kong.skp1.fa

cat OS.skp1.fa OBa.skp1.fa SB.skp1.fa > skp1.fa
perl ../DrawTree.pl --protein skp1.fa --align > log 2> log2 &
perl ../DrawTree.pl --protein kong.skp1.fa --align > log 2> log2 &

echo "get gene postion"
perl ../getidpos.pl --list PF01466.PF03931.tigr.pfam3.id --gff /home/jfchen/FFproject/seqlib/BGI_analysis_data/tigr_rice/tigr.all.final.gff --output OS.skp1.position

