perl matchlist.pl --gff ../input/RAP3.gff3.nr.gff --blasttable ../input/blasttable/OS_OS.blasttable --project rice
perl findTandemGeneDuplicates.pl -m 9 < rice.matchlist > rice.tandem.repeat.txt
perl formattandem.pl --tandem rice.tandem.repeat.txt --project rice > log 2> log2 &


