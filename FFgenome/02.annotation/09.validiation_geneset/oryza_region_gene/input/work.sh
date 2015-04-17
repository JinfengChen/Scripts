blastall -p blastn -i Oryza.fa -d OBa.glean.cds -e 1e-5 -o Oryza2glean.blastm8 -m 8 > log 2> log2 &

perl ../bin/OryzaName.pl --fasta MOC1_FJ032639.txt --region Moc1
perl ../bin/OryzaName.pl --fasta ADH1_FJ266021.txt --region Adh1
perl ../bin/OryzaName.pl --fasta Hd1_GQ407107.txt --region Hd1
cat *.final.fa > Oryza.final.fa
perl /home/jfchen/FFproject/tools/bin/runblastall_multicpu.pl -p blastn -i Oryza.final.fa -d OBa.gramene.FINAL_fgenesh.bestmodel.fa > log &
perl /home/jfchen/FFproject/tools/bin/runblastall_multicpu.pl -p blastn -i Oryza.final.fa -d OBa.glean.cds > log &
perl /home/jfchen/FFproject/tools/bin/runblastall_multicpu.pl -p blastn -i Oryza.final.fa -d Oryza_brachyantha.genome.super_scaffold.v1.0.fa > log &
perl /home/jfchen/FFproject/tools/bin/runblastall_multicpu.pl -p blastn -i Oryza.final.fa -d OBa.gramene.FINAL8_3_last.bestmodel.filter.fa > log 2> log2 &
perl /home/jfchen/FFproject/tools/bin/runblastall_multicpu.pl -p blastn -i Oryza.final.fa -d OBa.gramene.cds.fa > log 2> log2 &



