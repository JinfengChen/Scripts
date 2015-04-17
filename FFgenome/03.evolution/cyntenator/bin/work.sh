perl convertGFF.pl --gff ../input/tigr.all.final.gff --project Os
perl convertGFF.pl --gff ../input/OBa.all.chr.gff --project Ob
perl convertGFF.pl --gff ../input/Sb.final.gff --project Sb
perl convertBLAST4m8.pl --blast ../input/all_vs_all.blast.m8 --project 3way > log 2> log2 &
cd input
../bin/cyntenator -t "((Os.txt Ob.txt) Sb.txt)" -h blast 3way.blast  > log 2> log2 &


perl convertGFFchr.pl --gff ../input/tigr.all.final.gff --project Os > log 2> log2 &
perl convertGFFchr.pl --gff ../input/OBa.all.chr.gff --project Ob > log 2> log2 &
perl convertGFFchr.pl --gff ../input/Sb.final.gff --project Sb > log 2> log2 &


