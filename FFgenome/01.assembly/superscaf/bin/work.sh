perl superscaf.pl -fasta test.fa -link ../input/joinScaffold.txt -cut ../input/split.txt -project super > log 
perl superscaf.pl -fasta /share/raid12/chenjinfeng/FFgenome/genome/scaffold/20k/super-scaffold.fa -link ../input/joinScaffold.txt -cut ../input/split.txt -project OBa > log

perl superscafV2.pl -fasta /share/raid12/chenjinfeng/FFgenome/genome/scaffold/20k/super-scaffold.fa -link ../input/joinScaffold.txt -gff /share/raid12/chenjinfeng/seqlib/BGI_analysis_data/brachyantha/Oryza_brachyantha.genome.super_scaffold.v1.0.glean.gff.150_filter.2K.gff -cut ../input/split.txt -project OBaV2 > OBaV2.log 2> OBaV2.log2 &


echo "link gramene gff"
echo "final2"
perl superscafV2.pl --fasta ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --link ../input/joinScaffold.txt --gff ../input/OBa.gramene.gff --cut ../input/split.txt --project Gramene
cat Gramene.super.gff Gramene.super.unmaped.gff > Gramene.all.gff
echo "final6"
perl superscafV2.pl --fasta ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --link ../input/joinScaffold.txt --gff ../input/OBa.gramene.gff --cut ../input/split.txt --project Gramene
cat Gramene.super.gff Gramene.super.unmaped.gff > Gramene.all.gff

echo "final6_best2"
perl superscafV2.pl --fasta ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --link ../input/joinScaffold.txt --gff ../input/OBa.gramene.clean.best2.gff --cut ../input/split.txt --project Gramene

echo "final6_best3"
perl superscafV2.pl --fasta ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --link ../input/joinScaffold.txt --gff ../input/OBa.gramene.FINAL.bestmodel.gff --cut ../input/split.txt --project Gramene

echo "final8_3"
perl superscafV2.pl --fasta ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --link ../input/joinScaffold.txt --gff ../input/OBa.gramene.FINAL8_3.bestmodel.gff --cut ../input/split.txt --project Gramene.FINAL8_3

echo "final9"
perl superscafV2.pl --fasta ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --link ../input/joinScaffold.txt --gff ../input/OBa.gramene.FINAL9.bestmodel.gff --cut ../input/split.txt --project Gramene.FINAL9

echo "final8_2"
perl superscafV2.pl --fasta ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --link ../input/joinScaffold.txt --gff ../input/OBa.gramene.FINAL8_2.bestmodel.gff --cut ../input/split.txt --project Gramene.FINAL8_2

echo "final9a"
perl superscafV2.pl --fasta ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --link ../input/joinScaffold.txt --gff ../input/OBa.gramene.FINAL9a.bestmodel.gff --cut ../input/split.txt --project Gramene.FINAL9a

echo "final9_best"
perl superscafV2.pl --fasta ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --link ../input/joinScaffold.txt --gff ../input/OBa.gramene.FINAL9_best.bestmodel.gff --cut ../input/split.txt --project Gramene.FINAL9_best


echo "final8_3"
perl superscafV2.pl --fasta ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --link ../input/joinScaffold.txt --gff ../input/OBa.gramene.FINAL8_3_last.bestmodel.filter.gff --cut ../input/split.txt --project Gramene.FINAL8_3_last
perl superscafV3.pl --fasta ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --link ../input/joinScaffold.txt --gff ../input/OBa.gramene.FINAL8_3_last.bestmodel.filter.gff --cut ../input/split.txt --project Gramene.FINAL8_3_last

echo "assign chromosome code"
perl AssinGeneCode4chr.pl --gff Gramene.FINAL8_3_last.chr.gff > codetable
perl AssinGeneCode4scaffold.pl --gff Gramene.FINAL8_3_last.super.gff --code codetable
perl AssinGeneCode4scaffold.pl --gff ../input/OBa.gramene.FINAL8_3_last.bestmodel.filter.gff --code codetable

echo "produce a chr.scaffold table that can be used in drawFeature"
perl scaffold2chrTable.pl --fasta Gramene.FINAL8_3.super.fa --link ../input/joinScaffold.txt

echo "use gene colinarity and bes to correct some join and split"
perl superscafV3.pl -fasta /home/jfchen/FFproject/seqlib/BGI_analysis_data/brachyantha/superscaffold/Oryza_brachyantha.genome.super_scaffold.v1.0.fa -link ../input/joinScaffold.txt -gff /home/jfchen/FFproject/seqlib/BGI_analysis_data/brachyantha/superscaffold/Oryza_brachyantha.genome.super_scaffold.v1.0.glean.gff.150_filter.2K.gff -cut ../input/split.txt -project OBav3 > OBav3.log 2> OBav3.log2 &

echo "to raju"
perl superscafV3.pl --fasta ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --link ../input/joinScaffold.txt --gff ../input/OBa.gramene.FINAL8_3_last.bestmodel.filter.gff --cut ../input/split.txt --project Gramene > log 2> log2 &
perl AssinGeneCode4chr.pl --gff Gramene.chr.gff --project Gramene.v1 > codetable
perl AssinGeneCode4scaffold.pl --gff ../input/OBa.gramene.FINAL8_3_last.bestmodel.filter.gff --code codetable --project Gramene.v1 > log 2> log2 &

echo "glaberrima"
perl superscaf4glaberrima.pl --fasta ../input/glaberrima/glaberrima.fasta --link ../input/glaberrima/JoinScaffold_glaberrima.txt --cut ../input/split.txt --gff ../input/glaberrima/glaberrima.final.gff --project glaberrima > log 2> log2 &

echo "OBav3"
perl superscafV3.pl --fasta ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --link ../input/joinScaffold.txt --gff /home/jfchen/FFproject/seqlib/BGI_analysis_data/brachyantha/superscaffold/Oryza_brachyantha.genome.super_scaffold.v1.0.glean.gff.150_filter.2K.gff --cut ../input/split.txt --project OBav3 > log 2> log2 &

echo "OBav4"
perl superscafV3.pl --fasta ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --link ../input/joinScaffold.txt --gff /home/jfchen/FFproject/seqlib/BGI_analysis_data/brachyantha/superscaffold/Oryza_brachyantha.genome.super_scaffold.v1.0.glean.gff.150_filter.2K.gff --cut ../input/split.txt --project OBav4 > log 2> log2 &

echo "Gramene v1.4, correct some inversion by synteny analysis"
perl superscafV3.pl --fasta ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --link ../input/joinScaffold.txt --gff ../input/OBa.gramene.FINAL8_3_last.bestmodel.filter.gff --cut ../input/split.txt --project Gramene > log 2> log2 &
perl AssinGeneCode4chr.pl --gff Gramene.chr.gff --project Gramene.v1.4 > codetable
perl AssinGeneCode4scaffold.pl --gff ../input/OBa.gramene.FINAL8_3_last.bestmodel.filter.gff --code codetable --project Gramene.v1.4 > log 2> log2 &
perl superscafV3.pl --fasta ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --link ../input/joinScaffold.txt --gff ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.manual.TE.gff --cut ../input/split.txt --project Gramene.manual.TE > log 2> log2 &

echo "Glean_v1.4"
perl superscafV3.pl --fasta ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --link ../input/joinScaffold.txt --gff /home/jfchen/FFproject/seqlib/BGI_analysis_data/brachyantha/superscaffold/Oryza_brachyantha.genome.super_scaffold.v1.0.glean.gff.150_filter.2K.gff --cut ../input/split.txt --project Glean.v1.4 > log 2> log2 &


echo "Fbox reannotation"
perl superscafV3.pl --fasta ../input/Oryza_brachyantha.genome.super_scaffold.v1.0.fa --link ../input/joinScaffold.txt --gff ../input/final_v2.obra.gff --cut ../input/split.txt --project Fbox


