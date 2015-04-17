echo "Reference"
echo "Gene mean min max"
cat all.con.gff3.clear.density | perl numberStat.pl -type mean:min:max
echo "DNA te mean min max"
cat all.con.RepeatMasker.out.gff.DNA.density | perl numberStat.pl -type mean:min:max
echo "LTR te mean min max"
cat all.con.RepeatMasker.out.gff.LTR.density | perl numberStat.pl -type mean:min:max
echo "Other te mean min max"
cat all.con.RepeatMasker.out.gff.Other.density | perl numberStat.pl -type mean:min:max
echo "query"
echo "Gene mean min max"
cat ass.scafSeq.gapfill3.glean.gff.chr.density | perl numberStat.pl -type mean:min:max
echo "DNA te mean min max"
cat allTE.gff.clear.chr.DNA.density | perl numberStat.pl -type mean:min:max
echo "LTR te mean min max"
cat allTE.gff.clear.chr.LTR.density | perl numberStat.pl -type mean:min:max
echo "Other te mean min max"
cat allTE.gff.clear.chr.Other.density | perl numberStat.pl -type mean:min:max
echo "Root mean min max"
cat root.soap.chr.density | perl numberStat.pl -type mean:min:max
echo "Shoot mean min max"
cat shoot.soap.chr.density | perl numberStat.pl -type mean:min:max

