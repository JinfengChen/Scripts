perl SumValid.pl --blasttable ../input/Oryza.final.fa.BGIglean.blasttable > glean.summary.txt
perl SumValid.pl --blasttable ../input/Oryza.final.fa.gramene.blasttable > gramene.summary.txt

perl genbank2gff.pl --gb ../input/MOC1_FJ032639.gb --project Moc1 >log
perl genbank2gff.pl --gb ../input/ADH1_FJ266021.gb --project Adh1 >log
perl genbank2gff.pl --gb ../input/Hd1_GQ407107.gb --project Hd1 >log

perl GFF2embl.pl -gff Moc1.gff -embl Moc1.embl -fasta Moc1.fasta
perl GFF2embl.pl -gff Adh1.gff -embl Adh1.embl -fasta Adh1.fasta
perl GFF2embl.pl -gff Hd1.gff -embl Hd1.embl -fasta Hd1.fasta


perl GetSubGFF.pl --qryhead Adh1_Scaffold000030_3212299_3434574 --qrygff ../input/OBa.all.gff --qryfasta ../input/OBa.all.fa
perl GetSubGFF.pl --qryhead Moc1_Scaffold000018_3823279_4020080 --qrygff ../input/OBa.all.gff --qryfasta ../input/OBa.all.fa 
perl GetSubGFF.pl --qryhead Hd1_Scaffold000015_6859960_6970148 --qrygff ../input/OBa.all.gff --qryfasta ../input/OBa.all.fa 

perl runblast2seq.pl
perl run2act.pl

perl genbank2gff.pl --gb ../input/FQ378032.gb --project FQ378032
perl genbank2gff.pl --gb ../input/FQ378033.gb --project FQ378033
perl genbank2gff.pl --gb ../input/DQ810282.gb --project DQ810282
