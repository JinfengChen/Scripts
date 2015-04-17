perl GetPfamGeneID.pl --iprscan ../input/iprscan/OS.iprscan --domain PF00646
perl getidseq.pl -list PF00646.OS.id -fasta ../input/tigr.all.pep.final.fa -output tigr.fbox.fa > log 2> log2 &

phyml -i proteic -d aa > log 2> log2
perl /home/jfchen/software/tree/tree_plot.pl --png proteic_phyml_tree.txt

perl DrawTree.pl --protain protein.fa > log 2> log2 &


echo "get tigr fbox representive pep sequence"
perl getidseq.pl --list ./Fbox/FboxRice --fasta ../input/tigr.all.pep.final.fa --output tigr.fbox.representive.fa
echo "get PF00640 Fbox gene id for FF genome"
perl GetPfamGeneID.pl --iprscan ../input/iprscan/OB.iprscan --domain PF00646
echo "get pep sequence for FF fbox gene"
perl getidseq.pl --list PF00646.OB.id --fasta ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.glean.gff.150_filter.2K.pep --out OBa.fbox.fa
echo "create file for tree"
cat tigr.fbox.representive.fa OBa.fbox.fa > OBa.fbox4tree.fa
echo "align and draw tree"
perl DrawTree.pl --protein OBa.fbox4tree.fa --align > log 2> log2 &


echo "get gene contains two domain"
perl Get2PfamGeneID.pl --iprscan ../input/iprscan/OS.iprscan --domain1 PF00560 --domain2 PF00069
echo "get gene domains"
perl GetGenePfam.pl --iprscan ../input/iprscan/OS.iprscan --domain PF00646 > log 2> log2 &
#### iprscan miss some fbox gene, so we use hmmer3 hmmsearch to find domain
echo "search pfam with hmmer3 hmmsearch"
hmmsearch --tblout test.hmmer3.1 ../input/Pfam/Pfam-A.hmm tigr.fbox.representive.fa > log 2> log2 &
echo "get gene domain from hmmer hmmsearch results"
perl GetGenePfam4hmmer3.pl --hmmer3 test.hmmer3 --domain PF00646

echo "search pfam with hmmer3"
hmmsearch --tblout tigr.pfam3.hmmer3 ../input/Pfam/Pfam-A.hmm ../input/tigr.all.pep.final.fa > log 2> log2 &
hmmsearch --tblout OBa.pfam3.hmmer3 ../input/Pfam/Pfam-A.hmm ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.glean.gff.150_filter.2K.pep > log 2> log2 &
echo "get gene domain from hmmer results"
perl GetGenePfam4hmmer3.pl --hmmer3 tigr.pfam3.hmmer3 --domain PF00646 > log 2> log2 &
perl GetGenePfam4hmmer3.pl --hmmer3 OBa.pfam3.hmmer3 --domain PF00646 > log 2> log2 &
echo "get protein sequence of Fbox types"
perl GetTypeFasta.pl --type tigr.fbox.type --fasta ../input/tigr.all.pep.final.fa
perl GetTypeFasta.pl --type OBa.fbox.type --fasta ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.glean.gff.150_filter.2K.pep
echo "get protein sequence of ob and sb by hcluster"
perl GetTypeFastaRelated.pl --type tigr.fbox.type --fasta ../input/tigr.all.pep.final.fa --hcluster ../input/all_vs_all.blast.m8.solar.forHC.hcluster --ob ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.glean.gff.150_filter.2K.pep --sb ../input/Sb.final.pep > log1 2> log3 &
echo "add os"
perl GetTypeFastaRelated.pl --type tigr.fbox.type --fasta ../input/tigr.all.pep.final.fa --hcluster ../input/all_vs_all.blast.m8.solar.forHC.hcluster --ob ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.glean.gff.150_filter.2K.pep --sb ../input/Sb.final.pep --os ../input/tigr.all.pep.final.fa > log1 2> log3 &

echo "get related f box gene by type of three species and merge them together"
perl GetTypeFastaRelated4merge.pl --type tigr.fbox.type --fasta ../input/tigr.all.pep.final.fa --hcluster ../input/all_vs_all.blast.m8.solar.forHC.hcluster --ob ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.glean.gff.150_filter.2K.pep --sb ../input/Sb.final.pep --os ../input/tigr.all.pep.final.fa > log 2> log2 &

echo "align using hmmalign"
hmmalign --trim --outformat PSIBLAST ../input/hmm/F-box.hmm OBa.fbox4tree.fa > align.log 2> align.log2 &
echo "edit alignment in windows and draw tree using phyml"

