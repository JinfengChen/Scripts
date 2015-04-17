perl plotMethylation2chr.pl --chrlen ../input/IRGSP.chrlen --bar ../input/CHG.bar --project CHG_rice > log 2> log2 &
perl plotMethylation2chr.pl --chrlen ../input/IRGSP.chrlen --bar ../input/rice_methylation_bar --project rice > log 2> log2 &
perl plotMethylation2chr.pl --chrlen ../input/IRGSP.chrlen --bar ../input/rice_methylation_bar -type all --project rice > log 2> log2 &
perl plotFeature2chr.pl --chrlen ../input/IRGSP.chrlen --bar ../input/rice_feature_bar --type all --project rice > log 2> log2 &
perl plotFeature2chr.pl --chrlen ../input/IRGSP.chrlen --bar ../input/rice_feature_bar --epi --type all --project rice > log 2> log2 &
perl plotChipseq2chr.pl --chrlen ../input/IRGSP.chrlen --bar ../input/rice_chipseq_bar --type H3K4 --project rice_H3K4 > log 2> log2 &
perl plotChipseq2chr.pl --chrlen ../input/IRGSP.chrlen --bar ../input/rice_RPKM_bar --type RPKM --project Rice_RPKM > log 2> log2 &


echo "OBa"
perl splitTE.pl --gff ../input/OBa.all.fa.RepeatMasker.out.gff.chr
perl gff2windowsbar4TE_type.pl --chrlen ../input/OBa.chrlen --bar ../input/OBa_feature_bar --gff ../input/OBa.all.fa.RepeatMasker.out.gff.chr.type > log 2> log2 &

echo "IRGSP"
perl splitTE.pl --gff ../input/IRGSP.build5.RepeatMasker.out.gff.chr

echo "CentO for FF"
perl plotChipseq2chr.pl --chrlen ../input/OBa.add100.chrlen --bar ../input/Cent4542OBa.sort.bed.bar --type CentO --project OBa_CentO > log 2> log2 &

echo "CentF for FF"
perl plotChipseq2chr.pl --chrlen ../input/OBa.add100.chrlen --pointtype rec --bar ../input/ffcent2OBa.bed.bar --type CentF --project OBa_CentF > log 2> log2 &


