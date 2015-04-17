Test run: to get a chr.scaffold file, this file is then used to prepare chr dir of gff or RNAseq file.
perl ../bin/drawFeature.pl --dir ../data/axt/prenet_net_OBa.all2IRGSP --refseq ../data/seq/IRGSP.build5 --queryseq ../data/seq/OBa.all.fa -type lastz > log 2> log2 &

First run: to draw feature figure, but density should be modified according the result.
perl ../bin/drawFeature.pl --dir ../data/axt/prenet_net_OBa.all2IRGSP/ --refseq ../data/seq/IRGSP.build5 --queryseq ../data/seq/OBa.all.fa --qrytegff ../data/gff/OBa.all.fa.RepeatMasker.out.gff.chr --qrygenegff ../data/gff/OBa.all.gff.chr --qryroot ../data/gff/Root.bed.chr --qryshoot ../data/gff/Shoot.bed.chr --reftegff ../data/gff/IRGSP.build5.RepeatMasker.out.gff.chr --refgenegff ../data/gff/RAP3.gff3.nr.gff.chr --type lastz  --heat > log 2> log3 &

perl /share/raid12/chenjinfeng/tools/bin/qsub-sge.pl --resource vf=0.5G OBa2IRGSP2.sh &

Second And Third run: same as first run, after modifying color bin start and end.

Forth run: Add shoot transcritome data for reference, color table are modified according to Heapmap pipe by BGI.

perl /share/raid12/chenjinfeng/tools/bin/qsub-sge.pl --resource vf=0.5G OBa2IRGSP3.sh &

Fivth run: RNAseq data are measured by bp% of expressed sequence as well as these for gene and TE.

perl /share/raid12/chenjinfeng/tools/bin/qsub-sge.pl --resource vf=0.5G OBa2IRGSP3.sh &

Sixth run: Root transcriptone was removed for FF, and the first run on 159 server.

Senveth run: Methylation bed was added to figure.


