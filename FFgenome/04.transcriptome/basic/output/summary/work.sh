echo "test modified script to check if got same result with previous findoverlap method"
perl RNAseqDistribution.pl -bed ../OBa_GFF_gffjunction_tophat2/accepted_hits.bed -gene OBa.all.gff --repeat OBa.all.fa.RepeatMasker.out.gff --project OBa_BGI_Sum > log 2> log2 &

echo "Run summary on Gramene gene and manual repeat gff"
## remember this is based on bigscaffold which we use fpc conneted by super-scaffold
perl RNAseqDistribution.pl --bed ../OBa_GFF_alljunction_tophat2/accepted_hits.bed --gene Gramene.all.gff --repeat OBa.all.manual.TE.gff --project OBaGrameneGene > log 2> log2 &

perl RNAseqGene.pl --bed ../OBa_GFF_alljunction_tophat2/accepted_hits.bed --gene OBa.all.gff --project OBa.BGI

awk '$5 > 0' sumgene.sh.21171.qsub/sumgene.exon.coverage.table | wc -l

