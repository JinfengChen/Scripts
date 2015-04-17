perl mcscan_interval.pl -refgff ../input/Os.gff -qrygff ../input/Ob.gff -align ../input/ob_os.aligns > log &
perl mcscan_interval.pl -refgff ../input/Os.gff -refpep ../input/Os.pep -qrygff ../input/Ob.gff -qrypep ../input/Ob.pep -align ob_os.alians > log &
perl mcscan_interval.pl -refgff ../input/Os.gff -refpep ../input/Os.pep -qrygff ../input/Ob.gff -qrypep ../input/Ob.pep -align ../input/ob_os.aligns > log &

perl mcscan_interval.pl -refgff ../input/Os.gff -refpep ../input/Os.pep -refgff3 ../input/Os.gff3 -reffasta ../input/Os.fasta -qrygff ../input/Ob.gff -qrypep ../input/Ob.pep -qrygff3 ../input/Ob.gff3 -qryfasta ../input/Ob.fasta -align ob_os.aligns > log 2> log2 &
perl mcscan_interval.pl -refgff ../input/Os.gff -refpep ../input/Os.pep -refgff3 ../input/Os.gff3 -reffasta ../input/Os.fasta -qrygff ../input/Ob.gff -qrypep ../input/Ob.pep -qrygff3 ../input/Ob.gff3 -qryfasta ../input/Ob.fasta -align ../input/ob_os.aligns > log 2> log2 &

perl mcscan_coverage.pl --refgff ../input/Os.gff --qrygff ../input/Ob.gff --align ../input/ob_os.aligns --chrlen ../input/ob_os.chr > log

cat ../input/ob_os.mcl.align | grep "GLEAN" | cut -f 2 | sort | uniq | wc -l

