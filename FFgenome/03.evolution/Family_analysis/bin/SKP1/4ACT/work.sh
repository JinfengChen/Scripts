blastall -p tblastn -i flankgene.fa -d /home/jfchen/FFproject/seqlib/BGI_analysis_data/glaberrima/Oglaberrima_v1.0/glaberrima.fasta -o skp12gla.blast -e 1e-5 -m 8
perl /home/jfchen/FFproject/tools/bin/fastaDeal.pl --get_id Oglab07_0143 /home/jfchen/FFproject/seqlib/BGI_analysis_data/glaberrima/Oglaberrima_v1.0/glaberrima.fasta > Oglab07_0143.fa

blastall -p tblastn -i flankgene.fa -d /home/seqlib/plant_genomes/rice_9311/ChrAll.SupScaf -e 1e-5 -o skp12indica -m 8
perl /home/jfchen/FFproject/tools/bin/fastaDeal.pl --get_id Supscaffold8 /home/seqlib/plant_genomes/rice_9311/ChrAll.SupScaf > Supscaffold8.fa
perl /home/jfchen/FFproject/tools/bin/fastaDeal.pl --sub 10051897-10117644 Supscaffold8.fa > skp1_indica.txt
