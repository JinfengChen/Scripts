ln -s /share/raid12/chenjinfeng/FFgenome/drawFeature/data/gff/OBa.all.fa.RepeatMasker.out.gff.chr/ ./
ln -s /share/raid12/chenjinfeng/FFgenome/drawFeature/data/gff/IRGSP.build5.RepeatMasker.out.gff.chr ./
ln -s /share/raid12/chenjinfeng/FFgenome/drawFeature/data/gff/RAP3.gff3.nr.gff.chr/ ./
ln -s /share/raid12/chenjinfeng/FFgenome/drawFeature/data/gff/OBa.all.gff.chr/ ./

perl /share/raid12/chenjinfeng/tools/bin/fastaDeal.pl --attr id:len /share/raid12/chenjinfeng/FFgenome/drawFeature/data/seq/IRGSP.build5 > ricelen.txt
perl /share/raid12/chenjinfeng/tools/bin/fastaDeal.pl --attr id:len /share/raid12/chenjinfeng/FFgenome/drawFeature/data/seq/ffseq.chr > fflen.txt
sed 's/chr/A/' fflen.txt > fflen
sed 's/chr/A/' ricelen.txt > ricelen

cp /share/raid12/chenjinfeng/FFgenome/chralign/dotplot/bin/ob_os.blast ./
cat ob_os.blast | cut -f 1 | uniq > ffpara.gene
cat ob_os.blast | cut -f 2 | uniq > ricepara.gene
