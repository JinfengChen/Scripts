perl GFF2BED.pl -gff ../input/chr01/FF.Gene.gff -feature mRNA -bed FF.mRNA.bed
perl GFF2BED.pl -gff ../input/chr01/Rice.Gene.gff -feature mRNA -bed Rice.mRNA.bed

perl Expression2TE.pl -bed FF.mRNA.bed -gff ../input/chr01/FF.COPIA.gff -rpkm ../input/FF.shoot.rpkm -windows 10000 -project A01_Gene_COPIA_10000 > log
perl Expression2TE.pl -bed FF.mRNA.bed -gff ../input/chr01/FF.GYPSY.gff -rpkm ../input/FF.shoot.rpkm -windows 10000 -project A01_Gene_GYPSY_10000 > log
perl Expression2TE.pl -bed FF.mRNA.bed -gff ../input/chr01/FF.RT.gff -rpkm ../input/FF.shoot.rpkm -windows 10000 -project A01_Gene_RT_10000 > log
perl Expression2TE.pl -bed FF.mRNA.bed -gff ../input/chr01/FF.DNA.gff -rpkm ../input/FF.shoot.rpkm -windows 10000 -project A01_Gene_DNA_10000 > log
perl TE2species.pl -ortholog OB-OS.orth -bed1 FF.mRNA.bed -gff1 ../input/chr01/FF.COPIA.gff -bed2 Rice.mRNA.bed -gff2 ../input/chr01/Rice.COPIA.gff -windows 50000 -project A01_Gene_COPIA_2 > log
perl TE2species.pl -ortholog OB-OS.orth -bed1 FF.mRNA.bed -gff1 ../input/chr01/FF.GYPSY.gff -bed2 Rice.mRNA.bed -gff2 ../input/chr01/Rice.GYPSY.gff -windows 50000 -project A01_Gene_GYPSY_2 > log
perl TE2species.pl -ortholog OB-OS.orth -bed1 FF.mRNA.bed -gff1 ../input/chr01/FF.MITE.gff -bed2 Rice.mRNA.bed -gff2 ../input/chr01/Rice.MITE.gff -windows 50000 -project A01_Gene_MITE_2 > log
perl TE2species.pl -ortholog OB-OS.orth -bed1 FF.mRNA.bed -gff1 ../input/chr01/FF.CACTA.gff -bed2 Rice.mRNA.bed -gff2 ../input/chr01/Rice.CACTA.gff -windows 50000 -project A01_Gene_CACTA_2 > log
perl TE2species.pl -ortholog OB-OS.orth -bed1 FF.mRNA.bed -gff1 ../input/chr01/FF.MUDR.gff -bed2 Rice.mRNA.bed -gff2 ../input/chr01/Rice.MUDR.gff -windows 50000 -project A01_Gene_MUDR_2 > log


perl CompareRPKM.pl -bed FF.Gene.bed -rpkm1 ../input/FF.shoot.rpkm.bgi -rpkm2 ../input/FF.shoot.rpkm.tophat -project FF.compare

