perl rnaseq-search.pl --format bed ../input/OBa.gramene.FINAL8_2.bestmodel.gff ../input/accepted_hits.bed > log 2> log2 &
perl rnaseq-search.pl --format bowtie ../input/Gramene.chr.gff ../input/Gramene.RNAseq.bowtie > log 2> log2 &
