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
cat ass.scafSeq.gapfill3.RepeatMasker.denovo.change.out.gff.chr.DNA.density | perl numberStat.pl -type mean:min:max
echo "LTR te mean min max"
cat ass.scafSeq.gapfill3.RepeatMasker.denovo.change.out.gff.chr.LTR.density | perl numberStat.pl -type mean:min:max
echo "Other te mean min max"
cat ass.scafSeq.gapfill3.RepeatMasker.denovo.change.out.gff.chr.Other.density | perl numberStat.pl -type mean:min:max


