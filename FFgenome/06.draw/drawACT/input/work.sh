echo "split align"
perl ../bin/splitDAGalign.pl ob_os.matchlist.aligncoords > log 2> log2 &

echo "Tigr6"
perl ../bin/convertGFF.pl --gff tigr.all.final.gff --project Os
perl ../bin/splitGFF.pl tigr.all.final.gff 
perl ../bin/splitGFF.pl all.con.RepeatMasker.out.gff
perl ../bin/splitFASTA.pl all.con


echo "OBa"
perl ../bin/convertGFF.pl --gff OBa.chr.gff --project Ob
perl ../bin/splitGFF.pl OBa.chr.gff
perl ../bin/splitGFF.pl OBa.all.manual.TE.chr.gff
perl ../bin/splitFASTA.pl OBa.chr.fa

echo "gramenev1.4"
perl ../bin/splitGFF.pl Gramene.chr.gff
perl ../bin/splitGFF.pl Gramene.chr.manual.TE.gff
perl ../bin/splitFASTA.pl Gramene.chr.fa

echo "glaberrima"
perl ../bin/splitGFF.pl glaberrima.chr.gff 
perl ../bin/splitGFF.pl glaberrima.chr.fa.RepeatMasker.out.gff 
perl ../bin/splitFASTA.pl glaberrima.chr.fa

echo "9311"
perl ../bin/splitGFF.pl 9311.glean.final.gff 
perl ../bin/splitGFF.pl 9311.genome.final.fa.RepeatMasker.out.gff 
perl ../bin/splitFASTA.pl 9311.genome.final.fa

echo "epi"
perl ../bin/splitGFF.pl tigr6.1.H3K4me3.bed
perl ../bin/splitGFF.pl tigr6.1.mC.bed
perl ../bin/splitGFF.pl tigr6.1.shoot.bed

perl ../bin/splitGFF.pl gramenev1.4.H3K4me3.bed
perl ../bin/splitGFF.pl gramenev1.4.mC.bed
perl ../bin/splitGFF.pl gramenev1.4.shoot.bed

