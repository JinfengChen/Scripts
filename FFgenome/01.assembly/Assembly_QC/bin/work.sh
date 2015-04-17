perl verify2scafV4.pl -q ../input/BAC/BAC.fasta -t ../input/Gramene.chr.fa -b ../input/BAC2Gramenev1.4.blastm8 -i 99 -l 200 > log 2> log2 &
perl verify2scafV4.pl -q ../input/BAC/OB_B_PM_v6.txt -t ../input/Gramene.chr.fa -b ../input/Cent2Gramenev1.4.blastm8 -i 99 -l 200 > log 2> log2 &
perl verify2scafV4.pl -q ../input/BAC/OB_B_PM_v6.txt -t ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa -b ../input/Cent2Scaffold.blastm8 -i 99 -l 200 > log 2> log2 &


perl verify2scafV5.pl --query ../input/BAC/BAC.fasta -target ../input/Gramene.chr.fa -blastm8 ../input/BAC2Gramenev1.4.blastm8 --qryrepeatgff ../input/BAC/repeatmask/BAC.fasta.RepeatMasker.out.gff  --refgenegff ../input/Gramene.chr.gff --refrepeatgff ../input/Gramene.chr.manual.TE.gff -i 99 -l 200 > log 2> log2 &


perl statGFF.pl --fasta ../input/BAC/BAC.fasta --genegff ../input/BAC.gramene.gff.chr --tegff ../input/BAC.fasta.RepeatMasker.out.gff.chr > log


perl SumUnalign.pl --gene ../input/BAC.gramene.gff --repeat ../input/BAC.fasta.RepeatMasker.out.gff --unalign unalign.table --project BAC > log 2> log2 &


perl verify2scafV5.pl --query ../input/BAC/BAC.fasta -target ../input/Gramene.chr.fa -blastm8 ../input/BAC2Gramenev1.4.blastm8 --qrygenegff ../input/BAC.gramene.gff --qryrepeatgff ../input/BAC/repeatmask/BAC.fasta.RepeatMasker.out.gff  --refgenegff ../input/Gramene.chr.gff --refrepeatgff ../input/Gramene.chr.manual.TE.gff -i 99 -l 200 > log 2> log2 &


echo "ch1112 SD"

perl verify2scafV5.pl --query ../input/chr1112/chr1112.fasta -target ../input/Gramene.chr.fa -blastm8 ../input/chr11122Gramenev1.4.blastm8 --qrygenegff ../input/chr1112/chr1112.gff --qryrepeatgff ../input/chr1112/chr1112.fasta.RepeatMasker.out.gff  --refgenegff ../input/Gramene.chr.gff --refrepeatgff ../input/Gramene.chr.manual.TE.gff -i 99 -l 200 > log 2> log2 &

echo "chr03s"

perl verify2scafV5.pl --query ../input/chr03s/brachyantha.chr03s.fa -target ../input/chr03s/brachyanthaI.chr03s.fa -blastm8 ../input/chr03s2Gramenev1.4.chr03I.blastm8 --qrygenegff ../input/chr03s/brachyantha.chr03s.gene.gff --qryrepeatgff ../input/chr03s/brachyantha.chr03s.TE.gff  --refgenegff ../input/chr03s/brachyanthaI.chr03s.gene.gff --refrepeatgff ../input/chr03s/brachyanthaI.chr03s.TE.gff -i 99 -l 200 > log 2> log2 &

perl verify2scafV5.pl --query ../input/chr03s/brachyantha.chr03s.fa -target ../input/Gramene.chr.fa -blastm8 ../input/chr03s2Gramenev1.4.blastm8 --qrygenegff ../input/chr03s/brachyantha.chr03s.gene.gff --qryrepeatgff ../input/chr03s/brachyantha.chr03s.TE.gff --refgenegff ../input/Gramene.chr.gff --refrepeatgff ../input/Gramene.chr.manual.TE.gff -i 99 -l 200 > log 2> log2 &



