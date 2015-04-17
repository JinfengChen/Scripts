## small than 100 bp gaps:
perl ../../bin/gapfrequency.pl ../../data/scaffold20k.fa > gap.len

perl /share/raid12/chenjinfeng/tools/bin/fastaDeal.pl --attr id:gc gapflanking.fa > gapflanking.gc

perl ../../../../gcdraw/bin/slidingGC.pl /share/raid12/chenjinfeng/FFgenome/genome/scaffold/20k/super-scaffold.fa
mv /share/raid12/chenjinfeng/FFgenome/genome/scaffold/20k/super-scaffold.fa.gc  ./
gap2genome.r
gap20k.r

## samll than 200 bo gaps:

perl ../../bin/gapfrequency.pl ../../data/scaffold20k.fa > gap.len

perl /share/raid12/chenjinfeng/tools/bin/fastaDeal.pl --attr id:gc gapflanking.fa > gapflanking.gc

gap2genome.r
gap20k.r


