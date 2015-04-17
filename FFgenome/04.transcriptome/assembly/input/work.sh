blat Oryza_brachyantha.genome.super_scaffold.v1.0.fa rice.100scaftig.500.fa RNAseq2genome.500.blat.psl -noHead -minIdentity=95 > log 2> log2 &
blat Gramene.cds.fa rice.100scaftig.500.fa RNAseq2gene.500.blat.psl -noHead -minIdentity=95 > log 2> log2 &

cut -f10 RNAseq2gene.500.blat.psl | uniq | sort | uniq | wc -l
perl /home/jfchen/FFproject/tools/bin/bestAlign.pl RNAseq2gene.500.blat.psl --cutoff 0.3 > RNAseq2gene.500.blat.0.3.psl 
cut -f10 RNAseq2gene.500.blat.0.3.psl | uniq | sort | uniq | wc -l
cut -f10 RNAseq2genome.500.blat.psl | uniq | sort | uniq | wc -l
perl /home/jfchen/FFproject/tools/bin/bestAlign.pl RNAseq2genome.500.blat.psl --cutoff 0.3 > RNAseq2genome.500.blat.0.3.psl 
cut -f10 RNAseq2genome.500.blat.0.3.psl | uniq | sort | uniq | wc -l
perl /home/jfchen/FFproject/tools/bin/bestAlign.pl RNAseq2gene.500.blat.psl --cutoff 0.2 > RNAseq2gene.500.blat.0.2.psl 
cut -f10 RNAseq2gene.500.blat.0.2.psl | uniq | sort | uniq | wc -l

perl CoverageByTranscript.pl --transcript rice.100scaftig.500.fa  --gene Gramene.cds.fa --genome Oryza_brachyantha.genome.super_scaffold.v1.0.fa --project all > log 2> log2 &
perl selectBytranlaste.pl --fasta rice.100scaftig.500.fa --outfile rice.100scaftig.500.100aa.fa
perl CoverageByTranscript.pl --transcript rice.100scaftig.500.100aa.fa  --gene Gramene.cds.fa --genome Oryza_brachyantha.genome.super_scaffold.v1.0.fa --project 100aa > log 2> log2 &

perl /home/jfchen/FFproject/tools/bin/runblastall_multicpu.pl -p blastx -i missgene.fa -d /home/jfchen/FFproject/seqlib/BGI_analysis_data/plant_protein/plant_protein.fa > log 2> log2 &
perl /home/jfchen/FFproject/tools/bin/runRepeatMasker160.pl missgene.fa /home/jfchen/FFproject/FFgenome/02.annotation/01.repeat/1.repeat/manual/TElib-FF > log 2> log2 &

perl selectBymask.pl --masked rice.100scaftig.500.fa.RepeatMasker.masked --fasta rice.100scaftig.500.fa --outfile rice.100scaftig.500.nonrepeat.fa



