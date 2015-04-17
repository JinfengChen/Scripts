perl /home/jfchen/FFproject/tools/bin/getGene.pl --type exon ../input/Gramene.chr.gff ../input/Gramene.chr.fa >log 2>log2

blat ../input/all.con test.exon.fa test.psl > log 2> log2 &
perl /home/jfchen/FFproject/tools/bin/bestAlign.pl exon.psl --cutoff 0.5 > exon.best.psl
perl splitexon.pl --psl exon.best.psl > log 2> log2 &

echo "FF2rice"
perl /home/jfchen/FFproject/tools/bin/getGene.pl --type exon ../input/Gramene.chr.gff ../input/Gramene.chr.fa >log 2>log2

perl runblastall_multicpu_m8.pl -i Gramenev1.4.exon.fa -d ../input/all.con > log 2> log2 &

perl splitexon.pl --blastm8 Gramenev1.4.exon.fa.blast.m8 --repeat ../input/all.con.RepeatMasker.out.gff > log 2> log2 &

echo "rice2FF"
perl /home/jfchen/FFproject/tools/bin/getGene.pl --type exon ../input/tigr.all.final.gff ../input/all.con > tigr.exon.fa 2> log2 &

perl runblastall_multicpu_m8.pl -i tigr.exon.fa -d ../input/Gramene.chr.fa > log 2> log2 &

perl splitexon.pl --blastm8 tigr.exon.fa.blast.m8 --repeat ../input/Gramene.chr.manual.TE.gff > log 2> log2 &


echo "id"
perl getidseq.pl -l rice.id -f /home/jfchen/FFproject/seqlib/BGI_analysis_data/tigr_rice/all.pep -o rice.id.fa
blastall -p blastp -i rice.id.fa -d /home/seqlib/nr/nr -o rice2nr.blast.m8 -m 8 -e 1e-5



