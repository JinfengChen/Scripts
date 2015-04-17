echo "Reference"
echo "Gene mean min max"
cat RAP3.gff3.nr.gff.chr.density | perl numberStat.pl -type mean:min:max
echo "DNA te mean min max"
cat IRGSP.build5.RepeatMasker.out.gff.chr.DNA.density | perl numberStat.pl -type mean:min:max
echo "LTR te mean min max"
cat IRGSP.build5.RepeatMasker.out.gff.chr.LTR.density | perl numberStat.pl -type mean:min:max
echo "Other te mean min max"
cat IRGSP.build5.RepeatMasker.out.gff.chr.Other.density | perl numberStat.pl -type mean:min:max
echo "query"
echo "Gene mean min max"
cat OBa.all.gff.chr.density | perl numberStat.pl -type mean:min:max
echo "DNA te mean min max"
cat OBa.all.fa.RepeatMasker.out.gff.chr.DNA.density | perl numberStat.pl -type mean:min:max
echo "LTR te mean min max"
cat OBa.all.fa.RepeatMasker.out.gff.chr.LTR.density | perl numberStat.pl -type mean:min:max
echo "Other te mean min max"
cat OBa.all.fa.RepeatMasker.out.gff.chr.Other.density | perl numberStat.pl -type mean:min:max
echo "Root mean min max"
cat Root.bed.chr.density | perl numberStat.pl -type mean:min:max
echo "Shoot mean min max"
cat Root.bed.chr.density | perl numberStat.pl -type mean:min:max

