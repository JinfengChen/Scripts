### FF
###prepare gene
perl convertGFF2pos.pl -gff /home/jfchen/FFproject/seqlib/BGI_analysis_data/brachyantha/bigscaffold/OBa.super.gff -project FF_gene_bed > log 2> log2 &
###prepare CG,CHG,CHG bed file
perl count2bed_type.pl -coutdir ../input/FF_BS_methylation -type CG -bed FF_methylation_bed > log 2> log2 &
perl count2bed_type.pl -coutdir ../input/FF_BS_methylation -type CHG -bed FF_methylation_bed> log 2> log2 &
perl count2bed_type.pl -coutdir ../input/FF_BS_methylation -type CHH -bed FF_methylation_bed> log 2> log2 &
### calculate methylation level
perl methylevel.pl -gene FF_gene_bed -type CG -methy FF_methylation_bed --project FF > log 2> log2 &
perl methylevel.pl -gene FF_gene_bed -type CHG -methy FF_methylation_bed --project FF > log 2> log2 &
perl methylevel.pl -gene FF_gene_bed -type CHH -methy FF_methylation_bed --project FF > log 2> log2 &
### draw level figure
perl drawmethy.pl -CG FF.CG.level.sum -CHG FF.CHG.level.sum -CHH FF.CHH.level.sum -project FF
### combined the up two step into one
perl drawlevel.pl -bed FF_gene_bed -methy FF_methylation_bed -gene -project FF > log 2> log2 &
### prepare repeat
perl repeatGFF2pos.pl -gff /home/jfchen/FFproject/seqlib/BGI_analysis_data/brachyantha/bigscaffold/OBa.super.repeat.gff -project FF_repeat_bed > log 2> log2 &
### draw repeat level 
perl drawlevel.pl -bed FF_repeat_bed -methy FF_methylation_bed -repeat -project FF > log 2> log2 &

### draw qunitiles
perl drawqunitiles.pl -data ./FF -group ../input/Expression/FF.gene.group -type CG -project FF
### draw mCvsRPKM
perl bodymCvsRPKM.pl -rpkm ../input/Expression/FF.rpkm.tophat -body ./FF/FF.gene.CG.level.part -project FF


## rice
###prepare gene
perl convertGFF2pos.pl -gff /home/jfchen/FFproject/seqlib/BGI_analysis_data/rice_data/RAP3.gff3.nr.gff -project rice_gene_bed
###prepare CG,CHG,CHG bed file
perl count2bed_type.pl -coutdir ../input/IRGSP_BS_methylation -type CG -bed rice_methylation_bed > log 2> log2 &
perl count2bed_type.pl -coutdir ../input/IRGSP_BS_methylation -type CHG -bed rice_methylation_bed > log 2> log2 &
perl count2bed_type.pl -coutdir ../input/IRGSP_BS_methylation -type CHH -bed rice_methylation_bed > log 2> log2 &

### calculate methylation level
perl methylevel.pl -gene rice_gene_bed -type CG -methy rice_methylation_bed --project rice > log 2> log2 &
perl methylevel.pl -gene rice_gene_bed -type CHG -methy rice_methylation_bed --project rice > log 2> log2 &
perl methylevel.pl -gene rice_gene_bed -type CHH -methy rice_methylation_bed --project rice > log 2> log2 &
### draw level figure
perl drawmethy.pl -CG rice.CG.level.sum -CHG rice.CHG.level.sum -CHH rice.CHH.level.sum -project rice

### combined the up two step into one 
perl drawlevel.pl -bed rice_gene_bed -methy rice_methylation_bed -gene -project rice > log 2> log2 &

### prepare repeat 
perl repeatGFF2pos.pl -gff /home/jfchen/FFproject/FFgenome/06.draw/drawFeature/data/gff/IRGSP.build5.RepeatMasker.out.gff -project rice_repeat_bed > log 2> log2 &
### draw repeat level
perl drawlevel.pl -bed rice_repeat_bed -methy rice_methylation_bed -repeat -project rice > log 2> log2 &
perl drawlevelDZ.pl -bed rice_repeat_bed -methy rice_methylation_bed -gene -repeat -project rice > log 2> log2 &

###OBa 
###
### prepare repeat, manual repeat by dongying
perl repeatGFF2pos.pl --gff /home/jfchen/FFproject/seqlib/BGI_analysis_data/brachyantha/bigscaffold/OBa.all.manual.TE.chr.gff --project OBa_repeat_bed > log 2> log2 &
### draw repeat level


