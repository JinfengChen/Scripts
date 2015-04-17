perl /home/jfchen/FFproject/tools/bin/getGene.pl sub1.fregment.gff sub1.fregment.fa

perl /home/jfchen/FFproject/tools/bin/fastaDeal.pl --attr id sub1.fa > sub1.id
perl ../getidseq.pl --list sub1.OS.id --fasta /home/jfchen/FFproject/seqlib/BGI_analysis_data/tigr_rice/tigr.all.cds -o sub1.OS.cds.fa
perl ../getidseq.pl --list sub1.OBa.id --fasta /home/jfchen/FFproject/seqlib/BGI_analysis_data/brachyantha/superscaffold/Oryza_brachyantha.genome.super_scaffold.v1.0.glean.gff.150_filter.2K.cds > sub1.OBa.cds.fa
cat sub1A.cds.fa sub1.OS.cds.fa sub1.OBa.cds.fa > sub1all.cds.fa
cp sub1all.cds.fa ../../../adaptive_selection_analysis/input/

blastall -p blastn -i sub1.OS.cds.fa -d /home/jfchen/FFproject/seqlib/BGI_analysis_data/glaberrima/Oglaberrima_v1.0/glaberrima.fasta -o sub12glaberrima -e 1e-5 -m 8
perl /home/jfchen/FFproject/tools/bin/fastaDeal.pl --get_id Oglab09_0051 /home/jfchen/FFproject/seqlib/BGI_analysis_data/glaberrima/Oglaberrima_v1.0/Oglaberrima_chr9_v1.0.fasta > Oglab09_0051.fa

blastall -p blastn -i sub1.OS.cds.fa -d /home/seqlib/plant_genomes/rice_9311/ChrAll.SupScaf -e 1e-5 -o sub12indica -m 8

