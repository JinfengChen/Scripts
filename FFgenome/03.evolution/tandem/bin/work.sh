perl annotationTandem.pl --iprscan ../input/representative_orf.fa.nr.pep.iprscan --tandem OS.tandem_repeat_gene.txt --hcluster ../input/all_vs_all.blast.m8.solar.forHC.hcluster --project IRGSP > log 2> log2 &

echo "Tigr6 rice"
perl undirected.pl --blast ../input/blasttable/TR_TR.blasttable --gff ../input/tigr.all.final.gff > log 2> log2 &
perl annotationTandem.pl --iprscan ../input/tigr.all.final.iprscan --tandem TIGR6.tandem_repeat_gene.txt --hcluster ../input/tigr6.hcluster --ortholog ../input/TIGR6.ortholog.table --project TIGR6 > log 2> log2 &

