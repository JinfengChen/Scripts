perl inversion.pl --gff1 ../input/tigr.all.final.gff --gff2 ../input/Gramene.chr.gff --inversion ../input/inversion.v1.txt --blastm8 ../input/all_vs_all.blast.m8 | sort > summary.inversion

cat draw.r | R --vanilla --slave

perl invertrepeat.pl --inversion summary.inversion --refseq ../input/all.con --qryseq ../input/Gramene.chr.fa > log 2> log2 &

perl invertrepeat2.pl --inversion summary.inversion --refseq ../input/all.con --qryseq ../input/Gramene.chr.fa --refgff ../input/all.con.RepeatMasker.out.gff --qrygff ../input/Gramene.chr.manual.TE.gff > log 2> log2 &

perl summaryIR.pl --IR IR2.list


