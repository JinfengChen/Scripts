perl ../bin/sim4_to_gff3.pl --alignrate 0.5 --identity 0.9 map.sim4 > map.gff
perl compareGFF.pl --gff2 final_v2.obra.gff --gff1 map.gff > log 2> log2 &
perl finalGFF.pl --gff1 final_v2.obra.gff --gff2 map.gff > final.gff
perl /data/jfchen/Final_work/03.synteny/bin/checkID2.pl --id1 cds --id2 list2
msort -k 1,n4 final.gff > final.sort.gff




