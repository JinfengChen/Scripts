echo "split gene/TE gff into sub files by chromatin state"
perl splitChromatin.pl --chrlen ../input/OBa.chrlen --tegff ../input/OBa.all.fa.RepeatMasker.out.gff.chr --genegff ../input/OBa.all.gff.chr --segment ../input/OBa_methylation_wave --project OBa > OBa.log 2> OBa.log2 &
perl splitChromatin.pl --chrlen ../input/IRGSP.chrlen --tegff ../input/IRGSP.build5.RepeatMasker.out.gff.chr --genegff ../input/RAP3.gff3.nr.gff.chr --segment ../input/rice_methylation_wave --project rice > rice.log 2> rice.log2 &
perl splitChromatin.pl --chrlen ../input/OBa.chrlen --tegff ../input/OBa.all.manual.TE.gff.chr --genegff ../input/OBa.all.gff.chr --segment ../input/OBa_methylation_wave --project OBa_manual > OBa.log 2> OBa.log2 &

echo "prepare statistics for gene of euchromatin and heterchromatin to compare between species"
perl prepareChromatinInf.pl --chromatin OBa_manual >log 2>log2 &
perl prepareChromatinInf.pl --chromatin rice > log 2> log2 &

echo "compare and draw"
perl compareChromatin.pl --inf1 OBa_manual_chromatin --inf2 rice_chromatin --project OBa2rice > log 2> log2 &


echo "add rpkm,methylation and te density"
perl prepareChromatinInf.pl --chromatin OBa_manual --meth ../input/OBa.gene.CG.level.part --density ../input/OBa_50000.gene_rpkm_TE.density > log 2> log2 &
perl prepareChromatinInf.pl --chromatin rice --meth ../input/rice.gene.CG.level.part --density ../input/rice_50000.gene_rpkm_TE.density > log 2> log2 &

echo "compare and draw"
perl compareChromatin.pl --inf1 OBa_manual_chromatin --inf2 rice_chromatin --project OBa2rice > log 2> log2 &

