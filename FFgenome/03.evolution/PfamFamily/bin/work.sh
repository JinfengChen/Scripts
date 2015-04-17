perl getidseq.pl -l ../input/Pfam/PF00646_OB_219 -f /share/raid12/chenjinfeng/seqlib/BGI_analysis_data/brachyantha/Oryza_brachyantha.genome.super_scaffold.v1.0.glean.gff.150_filter.2K.pep -o PF00646_OB.fa
perl getidseq.pl -l ../input/Pfam/PF00646_TR_651 -f /share/raid12/chenjinfeng/seqlib/BGI_analysis_data/tigr_rice/tigr.all.pep.final.fa -o PF00646_TR.fa
perl getidseq.pl -l ../input/Pfam/PF00646_OS_433 -f /share/raid12/chenjinfeng/seqlib/BGI_analysis_data/rice_data/representative_orf.fa.nr.pep -o PF00646_OS.fa

echo "rice OBa";
echo "Get iprscan result for gene in ortholog/nonortholog table and summary data"
perl getgeneipr.pl --table ../input/nonorth/TIGR6.nonortholog.table --iprscan ../input/iprscan/OS.iprscan --project TIGR6.nonorth
perl getgeneipr.pl --table ../input/nonorth/TIGR6.syn.ortholog.table --iprscan ../input/iprscan/OS.iprscan --project TIGR6.syn.orth
perl getgeneipr.pl --table ../input/nonorth/TIGR6.nonsyn.ortholog.table --iprscan ../input/iprscan/OS.iprscan --project TIGR6.nonsyn.orth
perl SumPfam.pl --iprscan TIGR6
perl Size2synteny.pl --pfamsum TIGR6.Pfam.summary --project TIGR6

perl getgeneipr.pl --table ../input/nonorth/OBRACH.nonortholog.table --iprscan ../input/iprscan/OB.iprscan --project OBRACH.nonorth
perl getgeneipr.pl --table ../input/nonorth/OBRACH.syn.ortholog.table --iprscan ../input/iprscan/OB.iprscan --project OBRACH.syn.orth
perl getgeneipr.pl --table ../input/nonorth/OBRACH.nonsyn.ortholog.table --iprscan ../input/iprscan/OB.iprscan --project OBRACH.nonsyn.orth
perl SumPfam.pl --iprscan OBRACH
perl Size2synteny.pl --pfamsum OBRACH.Pfam.summary --project OBRACH

perl SumPfam.pl --iprscan TIGR6vsOBRACH/compare_nonorth
perl SumPfam.pl --iprscan TIGR6vsOBRACH/compare_synorth

echo "OBa sorghum"
perl getgeneipr.pl --table ../input/nonorth/OrthOrNot/OBRACH_SBICO/OBRACH.nonortholog.table --iprscan ../input/iprscan/OB.iprscan --project OBRACH.nonorth
perl getgeneipr.pl --table ../input/nonorth/SytnteyOrNot/OBRACH_SBICO/OBRACH.syn.ortholog.table --iprscan ../input/iprscan/OB.iprscan --project OBRACH.syn.orth
perl getgeneipr.pl --table ../input/nonorth/SytnteyOrNot/OBRACH_SBICO/OBRACH.nonsyn.ortholog.table --iprscan ../input/iprscan/OB.iprscan --project OBRACH.nonsyn.orth
perl SumPfam.pl --iprscan ./OBRACHvsSBCIO/OBRACH
mv OBRACHvsSBCIO/OBRACH.Pfam.summary OBRACHvsSBCIO/OBRACH/
perl Size2synteny.pl --pfamsum OBRACHvsSBCIO/OBRACH/OBRACH.Pfam.summary --project OBRACH

perl getgeneipr.pl --table ../input/nonorth/OrthOrNot/OBRACH_SBICO/SBICO.nonortholog.table --iprscan ../input/iprscan/SB.iprscan --project SBICO.nonorth
perl getgeneipr.pl --table ../input/nonorth/SytnteyOrNot/OBRACH_SBICO/SBICO.syn.ortholog.table --iprscan ../input/iprscan/SB.iprscan --project SBICO.syn.orth
perl getgeneipr.pl --table ../input/nonorth/SytnteyOrNot/OBRACH_SBICO/SBICO.nonsyn.ortholog.table --iprscan ../input/iprscan/SB.iprscan --project SBICO.nonsyn.orth
perl SumPfam.pl --iprscan ./OBRACHvsSBCIO/SBCIO
mv OBRACHvsSBCIO/SBCIO.Pfam.summary OBRACHvsSBCIO/SBCIO/
perl Size2synteny.pl --pfamsum OBRACHvsSBCIO/SBCIO/SBCIO.Pfam.summary --project SBCIO

echo "rice sorghum"
perl getgeneipr.pl --table ../input/nonorth/OrthOrNot/TIGR6_SBICO/TIGR6.nonortholog.table --iprscan ../input/iprscan/OS.iprscan --project TIGR6.nonorth
perl getgeneipr.pl --table ../input/nonorth/SytnteyOrNot/TIGR6_SBICO/TIGR6.syn.ortholog.table --iprscan ../input/iprscan/OS.iprscan --project TIGR6.syn.orth
perl getgeneipr.pl --table ../input/nonorth/SytnteyOrNot/TIGR6_SBICO/TIGR6.nonsyn.ortholog.table --iprscan ../input/iprscan/OS.iprscan --project TIGR6.nonsyn.orth
perl SumPfam.pl --iprscan ./TIGR6vsSBICO/TIGR6
mv ./TIGR6vsSBICO/TIGR6.Pfam.summary ./TIGR6vsSBICO/TIGR6/
perl Size2synteny.pl --pfamsum ./TIGR6vsSBICO/TIGR6/TIGR6.Pfam.summary --project TIGR6


perl getgeneipr.pl --table ../input/nonorth/SytnteyOrNot/TIGR6_SBICO/SBICO.nonsyn.ortholog.table --iprscan ../input/iprscan/SB.iprscan --project SBICO.nonsyn.orth
perl getgeneipr.pl --table ../input/nonorth/SytnteyOrNot/TIGR6_SBICO/SBICO.syn.ortholog.table --iprscan ../input/iprscan/SB.iprscan --project SBICO.syn.orth
perl getgeneipr.pl --table ../input/nonorth/OrthOrNot/TIGR6_SBICO/SBICO.nonortholog.table --iprscan ../input/iprscan/SB.iprscan --project SBICO.nonorth




