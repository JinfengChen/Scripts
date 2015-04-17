####rice.Gene.rpkm.xls is FF transcriptome result from BGI.
awk -F"\t" '{if($1~/OBR/){print $1,$8}}' rice.Gene.rpkm.xls > FF.shoot.rpkm.bgi
cp /share/raid12/chenjinfeng/FFgenome/transcriptome/output/OBaShoot_GFF_tophat2/OBa.all.tophat.gff.expr FF.shoot.rpkm.tophat
cp /share/raid12/chenjinfeng/FFgenome/transcriptome/output/IRGSP_Shoot_GFF_tophat2/RAP3.gff3.nr.tophat.gff.expr ./RAP3.shoot.rpkm.tophat


