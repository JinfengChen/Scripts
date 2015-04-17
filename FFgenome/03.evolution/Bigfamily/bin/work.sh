/share/raid12/chenjinfeng/tools/hmmer3.0/bin/hmmscan -o Fbox.hmmscan --tblout Fbox.gene --domtblout Fbox.domain ../input/hmm/F-box.hmm ../input/PF00646_OB.fa > log 2> log2 &
perl DomainInGenome.pl -i ../input/PF00646_OB_cds.fa -o cds.domain > log 2> log2 &
perl /share/raid12/chenjinfeng/tools/bin/qsub-sge.pl --resource vf=0.1G domain.sh &
cat cds.domain | cut -f 1 | uniq | wc -l

perl FindDomainGene.pl -gff /share/raid12/chenjinfeng/seqlib/BGI_analysis_data/brachyantha/Oryza_brachyantha.genome.super_scaffold.v1.0.glean.gff.150_filter.2K.gff -domain test.domain -iprscan /share/raid12/chenjinfeng/seqlib/BGI_analysis_data/brachyantha/Oryza_brachyantha.genome.super_scaffold.v1.0.glean.gff.150_filter.2K.pep.iprscan > log 2> log2 &


