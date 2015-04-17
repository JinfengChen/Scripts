echo "R gene contain NBS domain,NB-ARC/PF00931"
perl ../Get2PfamGeneID4hmmer3.pl --hmmer3 ../tigr.pfam3.hmmer3 --domain1 PF00931 --domain2 PF00069 > NBS.rice
perl ../GetGenePfam4hmmer3NBS.pl --hmmer3 ../tigr.pfam3.hmmer3 --domain PF00931 > tigr.type

perl ../Get2PfamGeneID4hmmer3.pl --hmmer3 ../OBa.pfam3.hmmer3 --domain1 PF00931 --domain2 PF00069 > NBS.OBa
perl ../GetGenePfam4hmmer3NBS.pl --hmmer3 ../OBa.pfam3.hmmer3 --domain PF00931 > OBa.type

perl ../Get2PfamGeneID4hmmer3.pl --hmmer3 ../Sb.pfam3.hmmer3 --domain1 PF00931 --domain2 PF00069 > NBS.SB
perl ../GetGenePfam4hmmer3NBS.pl --hmmer3 ../Sb.pfam3.hmmer3 --domain PF00931 > SB.type


perl ../getidseq.pl --list NBS.rice --fasta ../../input/tigr.all.pep.final.fa --output NBS.rice.fa
hmmsearch --tblout LRR.rice.hmmsearch1 ../../input/Pfam/LRR.hmm NBS.rice.fa > LRR.rice.hmmsearch.out1
cut -d " " -f1 LRR.rice.hmmsearch1 | sort | uniq | wc -l


/home/biosoftware/hmmer-2.3.2/bin/hmmsearch ../../input/Pfam/LRR.smart.hmm NBS.rice.fa > NBS.rice.hmmsearch2
hmmpfam --compat ../../input/Pfam/LRR.smart.hmm NBS.rice.fa > NBS.rice.hmmsearch2

perl hmmer2table.pl --hmmer NBS.rice.hmmsearch2 > log
cp /home/biosoftware/iprscan/data/new_coil.mat ./new.mat

echo "final pipeline"
perl ../Get2PfamGeneID4hmmer3.pl --hmmer3 ../tigr.pfam3.hmmer3 --domain1 PF00931 --domain2 PF00069 > NBS.rice
perl ../getidseq.pl --list NBS.rice --fasta ../../input/tigr.all.pep.final.fa --output NBS.rice.fa
perl classifyNBS.pl --fasta NBS.rice.fa --project rice
perl getidpos.pl --list rice.type --gff /home/jfchen/FFproject/seqlib/BGI_analysis_data/tigr_rice/tigr.all.final.gff --output NBS.rice.inf
perl classifyNBSv2.pl --fasta NBS.rice.fa --project tigrv2 > tigrv2.summary &

perl ../Get2PfamGeneID4hmmer3.pl --hmmer3 ../OBa.pfam3.hmmer3 --domain1 PF00931 --domain2 PF00069 > NBS.OBa
perl ../getidseq.pl --list NBS.OBa --fasta ../../input/Oryza_brachyantha.genome.super_scaffold.v1.0.glean.gff.150_filter.2K.pep --output NBS.OBa.fa
perl classifyNBS.pl --fasta NBS.OBa.fa --project OBa > OBa.summary

perl ../Get2PfamGeneID4hmmer3.pl --hmmer3 ../Sb.pfam3.hmmer3 --domain1 PF00931 --domain2 PF00069 > NBS.SB
perl ../getidseq.pl --list NBS.SB --fasta ../../input/Sb.final.pep --output NBS.SB.fa
perl classifyNBS.pl --fasta NBS.SB.fa --project SB > SB.summary



echo "align"
hmmalign --trim --outformat PSIBLAST ../../input/Pfam/NB-ARC.v2.hmm tigr6.NBS.fulldomain.fa > log 2> log2 &



