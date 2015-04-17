perl getGene.pl ../input/OBa_final.final5.gff ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa -type mrna > log 2> log2 &

echo "ensembl gff contain many duplicted gene because job_no_retry, so we remove these duplicated copy before analysis"
perl formatGFF.pl --gff ../input/OBa_final.final5.gff
#perl GFF2BED.pl --gff OBa.gramene.gff --feature mRNA -bed OBa.gramene.bed
#cut -f2 OBa.gramene.bed | sort | uniq | wc -l
perl getGene.pl OBa.gramene.gff ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --type mrna > log > log2 &
perl check.pl --gff OBa.gramene.gff --fasta OBa.gramene.cds > id

echo "link gff into big-scaffold to valid by RNAseq"
ln -s /home/jfchen/FFproject/FFgenome/01.assembly/superscaf/bin/Gramene.all.gff ./

echo "gramene_v1 final1"
perl formatGFF4Gramene.pl --gff ../input/OBa_gramene_v1_final1.gff > log 2> log2 &
perl getGene.pl OBa.gramene.gff ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --type mrna > OBa.gramene.cds

perl formatGFF4Gramene.pl --gff ../input/OBa_gramene_v1_final2.gff
perl getGene.pl OBa.gramene.gff ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --type mrna > OBa.gramene.cds

perl formatGFF4Gramene.pl --gff ../input/OBa_gramene_v1_final6.gff
perl getGene.pl OBa.gramene.gff ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --type mrna > OBa.gramene.cds

echo "get longest and all transcripts"
perl formatGFF4Gramene.pl --gff ../input/OBa_gramene_v1_final6.gff > log
perl getGene.pl OBa.gramene.gff ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --type mrna > OBa.gramene.cds
perl getGene.pl OBa.gramene.alltranscript.gff ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --type mrna > OBa.gramene.alltranscript.cds

