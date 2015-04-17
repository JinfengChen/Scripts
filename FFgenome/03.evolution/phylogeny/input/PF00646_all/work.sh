/share/raid1/genome/bin/clustalw PF00646.fa > log 2> log2 &
perl ../bin/AlignIO.pl -informat clustalw -outformat fasta -inaln PF00646.aln -outaln PF00646.fasta
perl /share/raid12/chenjinfeng/tools/fasta2/fasta2nex.pl
cat PF00646.nex ./treemethod/bays.protein > PF00646.bays.nex
perl /share/raid12/chenjinfeng/tools/bin/qsub-sge.pl --resource vf=0.2G bays.sh &
/share/raid12/chenjinfeng/tools/FastTree/FastTree PF00646.fasta > PF00646.tree 2> log2 &
