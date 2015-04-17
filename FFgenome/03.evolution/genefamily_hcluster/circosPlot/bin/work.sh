echo "Prepare for all ortholog pair"
perl prepareCircos.pl --ortholog ../input/OBRACH_TIGR6.distance.txt --gff1 ../input/OBa.chr.gff --gff2 ../input/tigr.all.final.gff --chrlen1 ../input/OBa.chrlen --chrlen2 ../input/tigr.chrlen > log 2> log2 &

echo "Split ortholog pair into synteny and nonsynteny"
perl SplitDistance.pl --syntid ../input/TIGR.syn.ortholog.table --distance ../input/OBRACH_TIGR6.distance.txt

echo "Prepare for synteny and nonsyteny seperately"
perl prepareCircos.pl --ortholog ../input/OBRACH_TIGR6.syn.distance.txt --gff1 ../input/OBa.chr.gff --gff2 ../input/tigr.all.final.gff --chrlen1 ../input/OBa.chrlen --chrlen2 ../input/tigr.chrlen > log 2> log2 &
perl prepareCircos.pl --ortholog ../input/OBRACH_TIGR6.nonsyn.distance.txt --gff1 ../input/OBa.chr.gff --gff2 ../input/tigr.all.final.gff --chrlen1 ../input/OBa.chrlen --chrlen2 ../input/tigr.chrlen > log 2> log2 &


echo "find out these synteny ortholog that have duplicated copy on other chromosome"
perl FindDupliOrt.pl --gff1 ../input/OBa.chr.gff --gff2 ../input/tigr.all.final.gff --orth ../input/OBRACH_TIGR6.syn.distance.txt --syntid ../input/TIGR.syn.ortholog.table > log 2> log2 &
                                                                             
