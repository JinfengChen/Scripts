perl checkgap.pl --fasta ../OBa.gramene.cds.fa > gapgene.id2
perl getidgff.pl --list gapgene.id -gff ../Gramene.chr.gff -o gapgene.gff
cut -f2 gapgene.id | perl /home/jfchen/FFproject/tools/bin/numberStat.pl
perl selectGene.pl --id gapgene.id2 > checklist
perl prepareGene.pl --list checklist --fasta ../OBa.gramene.pep.fa --output gapgene.pep.fa
blastall -p blastp -i gapgene.pep.fa -d /home/jfchen/FFproject/seqlib/BGI_analysis_data/plant_protein/plant_protein.fa -e 1e-5 -o gapgene2plantprotein.blastp &
perl /home/jfchen/FFproject/tools/bin/blast_parser.pl --tophit 1 gapgene2plantprotein.blastp > gapgene2plantprotein.blastp.table

blastall -p blastp -i gapgene.pep.fa -d /home/jfchen/FFproject/seqlib/BGI_analysis_data/tigr_rice/tigr.all.pep.final.fa -e 1e-5 -o gapgene2tigr.blastp &
perl /home/jfchen/FFproject/tools/bin/blast_parser.pl --tophit 1 -topmatch 1 gapgene2tigr.blastp > gapgene2tigr.blastp.table


perl gapgenehit.pl --pep gapgene.pep.fa --gff gapgene.gff --blast gapgene2tigr.blastp.table --list checklist





#version2, filter out opposite strand gene

perl checkgap.pl --fasta ../Gramene.cds.fa > gapgene.id
perl getidgff.pl --list gapgene.id -gff ../Gramene.chr.gff -o gapgene.gff
cut -f2 gapgene.id | perl /home/jfchen/FFproject/tools/bin/numberStat.pl > gapstat.txt
perl selectGene.pl --id gapgene.id > checklist
perl prepareGene.pl --list checklist --fasta ../Gramene.pep.fa --output gapgene.pep.fa
blastall -p blastp -i gapgene.pep.fa -d /home/jfchen/FFproject/seqlib/BGI_analysis_data/tigr_rice/tigr.all.pep.final.fa -e 1e-5 -o gapgene2tigr.blastp &


perl /home/jfchen/FFproject/tools/bin/blast_parser.pl --tophit 1 -topmatch 1 gapgene2tigr.blastp > gapgene2tigr.blastp.table
perl gapgenehit.pl --pep gapgene.pep.fa --gff gapgene.gff --blast gapgene2tigr.blastp.table --list checklist > summary.100

blastall -p blastp -i gapgene.pep.all.fa -d /home/jfchen/FFproject/seqlib/BGI_analysis_data/tigr_rice/tigr.all.pep.final.fa -e 1e-5 -o gapgeneall2tigr.blastp &


perl /home/jfchen/FFproject/tools/bin/blast_parser.pl --tophit 1 -topmatch 1 gapgeneall2tigr.blastp > gapgeneall2tigr.blastp.table
perl gapgenehit.pl --pep gapgene.pep.all.fa --gff gapgene.gff --blast gapgeneall2tigr.blastp.table --list checklist.all > summary.all

awk '$5 > 0' summary.100 > summary.100.nonzero
awk '$5 > 0' summary.all > summary.all.nonzero
cat coverage.r | R --vanilla --slave


