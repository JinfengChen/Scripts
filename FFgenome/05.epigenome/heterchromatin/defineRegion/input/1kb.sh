perl ../bin/bed2windowsbar4methylation_2.pl -window 10000 -step 10000 -chrlen IRGSP.chrlen -bar rice_feature_bar_10kb -bed rice_methylation_bed > log 2> log2
perl ../bin/gff2windowsbar4gene.pl -window 10000 -step 10000 -chrlen IRGSP.chrlen -bar rice_feature_bar_10kb -gff RAP3.gff3.nr.gff.chr > log 2> log2
perl ../bin/gff2windowsbar4TE.pl -window 10000 -step 10000 -chrlen IRGSP.chrlen -bar rice_feature_bar_10kb -gff IRGSP.build5.RepeatMasker.out.gff.chr > log 2> log2
