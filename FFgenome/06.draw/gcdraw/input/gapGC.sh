perl /share/raid12/chenjinfeng/tools/bin/getlargeseq.pl -i /share/raid12/chenjinfeng/FFgenome/genome/scaffold/20k/scaffold.fa -o largescaffold.fa -l 1000000
perl ../bin/slidingGC.pl largescaffold.fa > log
perl /share/raid12/chenjinfeng/FFgenome/assmbly/denovostat/bin/gapfrequency.pl largescaffold.fa > log 
mv gaps/*.gaps ./GCcontent/
cd GCcontent/
perl ../../bin/slidingGC.pl > log
