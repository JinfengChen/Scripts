perl pair_end_in_scaffold.pl ALL scaffold.lst ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa  soap.lst > log 2> log2 &

perl pair_end_in_scaffold.pl ALL scaffold.lst ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa  soap1.lst > log 2> log2 &

perl pair_end_in_scaffoldv2.pl ALL scaffold1.lst ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa soap1.lst > log 2> log2 &

echo "revision, chr11/12, chr04"
perl pair_end_in_scaffoldv2.pl PRE chr.lst ../input/Gramene.chr.fa soaprefchr.lst > log 2> log2 &
echo "modify scaff_config_pos.log and run DRAW"
perl pair_end_in_scaffoldv2.pl DRAW chr.lst ../input/Gramene.chr.fa soaprefchr.lst > log 2> log2 &


echo "get sub soap for draw"
perl getsubsoap.pl --soap soaprefchr.lst --list revision.mapping.pos --out ../input/soaprefdraw > log 2> log2 &
echo "generate draw scaffold list"
perl /home/jfchen/159/FFproject/tools/bin/fastaDeal.pl --attr id ../input/draw.fasta > draw.scaffold.lst
echo "generate draw soap list"
ls /home/jfchen/159/FFproject/FFgenome/01.assembly/PEscafold/input/soaprefdraw/*.soap > draw.soap.lst
echo "draw"
perl pair_end_in_scaffoldv2.pl ALL draw.scaffold.lst ../input/draw.fasta draw.soap.lst > log 2> log2 &
perl pair_end_in_scaffoldv2.pl DRAW draw.scaffold.lst ../input/draw.fasta draw.soap.lst > log 2> log2 &

echo "coverage"
/home/jfchen/159/revision/coverage/2.7.7/soap.coverage -cvg -il draw.soap.lst -refsingle ../input/draw.fasta -o result.txt -depthsingle all.coverage -plot distribution.txt 0 1000

/home/jfchen/159/revision/coverage/2.7.7/soap.coverage -cvg -il draw.soap.lst -refsingle ../input/draw.fasta -o result.txt -depthsingle all.coverage -window 10bpwin.coverage 10


echo "add coverage to PE figure"
perl pair_end_in_scaffoldv2check.pl DRAW draw.scaffold.lst ../input/draw.fasta draw.soap.lst all.coverage > log 2> log2 &

perl pair_end_in_scaffoldv2check.pl DRAW draw.scaffold.lst ../input/draw.fasta draw.soap.lst 10bpwin.coverage > log 2> log2 &

echo "run 4 chr11,12 sd"
perl pair_end_in_scaffoldv2check4chr11sd.pl DRAW draw.scaffold.lst ../input/draw.fasta draw.soap.lst 10bpwin.coverage > log 2> log2 &

perl pair_end_in_scaffoldv2check4chr12sd.pl DRAW draw.scaffold.lst ../input/draw.fasta draw.soap.lst 10bpwin.coverage > log 2> log2 &

echo "draw for inversion SV"
perl getsubsoap.pl --soap soaprefchr.lst --list inversion_SV.pos --out ../input/soaprefinversion > log 2> log2 &
perl /home/jfchen/159/FFproject/tools/bin/fastaDeal.pl --attr id ../input/inversionSV200kb.fasta > draw.inversionSV200kb.lst
ls /home/jfchen/159/FFproject/FFgenome/01.assembly/PEscafold/input/soaprefinversion/*.soap > draw.inversionSV200kb.soap.lst

/home/jfchen/159/revision/coverage/2.7.7/soap.coverage -cvg -il draw.inversionSV200kb.soap.lst -refsingle ../input/inversionSV200kb.fasta -o inversion_SV200kb.result.txt -depthsingle inversion_SV200kb.all.coverage -window inversion_SV200kb.100bpwin.coverage 100

perl pair_end_in_scaffoldv2check4SV.pl ALL draw.inversionSV200kb.scaffold.lst ../input/inversionSV200kb.fasta draw.inversionSV200kb.soap.lst inversion_SV200kb.100bpwin.coverage > log 2> log2 &
perl pair_end_in_scaffoldv2check4SV.pl DRAW draw.inversionSV200kb.scaffold.lst ../input/inversionSV200kb.fasta draw.inversionSV200kb.soap.lst inversion_SV200kb.100bpwin.coverage > log 6> log2 &

echo "H1 draw"
perl getsubsoap.pl --soap soaprefchr.lst --list H1.pos --out ../input/soaprefH1 > log 2> log2 &
perl /home/jfchen/159/FFproject/tools/bin/fastaDeal.pl --attr id ../input/H1.fasta > draw.H1.lst
ls /home/jfchen/159/FFproject/FFgenome/01.assembly/PEscafold/input/soaprefH1/*.soap > draw.H1.soap.lst

/home/jfchen/159/revision/coverage/2.7.7/soap.coverage -cvg -il draw.H1.soap.lst -refsingle ../input/H1.fasta -o H1.result.txt -depthsingle H1.all.coverage -window H1.100bpwin.coverage 100

perl pair_end_in_scaffoldv2check.pl ALL draw.H1.lst ../input/H1.fasta draw.H1.soap.lst H1.100bpwin.coverage > log 2> log2 &

echo "draw heterochromatin"
perl pair_end_in_scaffoldv2check4H1.pl ALL draw.H1.lst ../input/H1.fasta draw.H1.soap.lst H1.100bpwin.coverage > log 2> log2
cp OB_chr04_1057470_1668757.pdf /home/jfchen/159/Data/Final_work/09.dotplotRegion/bin/DrawN4heterochromatin/H1/

perl pair_end_in_scaffoldv2check4H7.pl ALL draw.H7H8.lst ../input/draw.fasta draw.soap.lst 10bpwin.coverage > log 2> log2 &
cp OB_chr04_6039164_6768670.* /home/jfchen/159/Data/Final_work/09.dotplotRegion/bin/DrawN4heterochromatin/H7/

perl pair_end_in_scaffoldv2check4H8.pl ALL draw.H7H8.lst ../input/draw.fasta draw.soap.lst 10bpwin.coverage > log 2> log2 &
cp OB_chr04_7027922_7396967.* /home/jfchen/159/Data/Final_work/09.dotplotRegion/bin/DrawN4heterochromatin/H8/


echo "chr03 short arm"
ls /home/jfchen/159/FFproject/FFgenome/01.assembly/PEscafold/input/soaprefchr03s/*.soap > draw.chr03s.soap.lst
perl pair_end_in_scaffoldv2check4chr03s.pl ALL draw.chr03s.lst ../input/brachyantha.chr03s.fa draw.chr03s.soap.lst > log 2> log2 &
perl pair_end_in_scaffoldv2check4chr03s.pl DRAW draw.chr03s.lst ../input/brachyantha.chr03s.fa draw.chr03s.soap.lst > log 2> log2 &
perl pair_end_in_scaffold.pl DRAW draw.chr03s.lst ../input/brachyantha.chr03s.fa draw.chr03s.soap.lst > log 2> log2 &

ls /home/jfchen/159/FFproject/FFgenome/01.assembly/PEscafold/input/soaprefchr03sII/*.soap > draw.chr03sII.soap.lst
perl pair_end_in_scaffold.pl ALL draw.chr03sII.lst ../input/brachyanthaII.chr03s.fa draw.chr03sII.soap.lst > log 2> log2 &

echo "chr03 short arm regions"
perl getsubsoap.pl --soap draw.chr03s.soap.lst --list chr03s_regions.pos --out ../input/soaprefchr03s_regions > log 2> log2 &
perl /home/jfchen/159/FFproject/tools/bin/fastaDeal.pl --attr id ../input/chr03s.fa > draw.chr03s.regions.lst
ls /home/jfchen/159/FFproject/FFgenome/01.assembly/PEscafold/input/soaprefchr03s_regions/*.soap > draw.chr03s.regions.soap.lst
perl pair_end_in_scaffoldv2check4chr03s2.5Mb.pl ALL draw.chr03s.regions.lst ../input/chr03s.fa draw.chr03s.regions.soap.lst > log 2> log2 &
perl pair_end_in_scaffoldv2check4chr03s3.4Mb.pl DRAW draw.chr03s.regions.lst ../input/chr03s.fa draw.chr03s.regions.soap.lst > log 2> log2 &

perl getsubsoap.pl --soap draw.chr03sII.soap.lst --list chr03sII_regions.pos --out ../input/soaprefchr03sII_regions > log 2> log2 &
perl /home/jfchen/159/FFproject/tools/bin/fastaDeal.pl --attr id ../input/chr03sII.fa > draw.chr03sII.regions.lst
ls /home/jfchen/159/FFproject/FFgenome/01.assembly/PEscafold/input/soaprefchr03sII_regions/*.soap > draw.chr03sII.regions.soap.lst
perl pair_end_in_scaffoldv2check4chr03s2.5Mb.pl ALL draw.chr03sII.regions.lst ../input/chr03sII.fa draw.chr03sII.regions.soap.lst > log 2> log2 &
perl pair_end_in_scaffoldv2check4chr03s3.4Mb.pl DRAW draw.chr03sII.regions.lst ../input/chr03sII.fa draw.chr03sII.regions.soap.lst > log 2> log2 &




