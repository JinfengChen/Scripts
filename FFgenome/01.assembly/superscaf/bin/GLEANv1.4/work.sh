perl gff2table.pl --gff Glean.v1.4.chr.gff > Glean.mRNA.chr.table
perl gff2table.pl --gff ../Gramenev1.4/Gramene.chr.gff > Gramene.mRNA.chr.table
perl /home/jfchen/FFproject/tools/findOverlap/findOverlap.pl Glean.mRNA.chr.table Gramene.mRNA.chr.table > Glean2Gramene

cat fbox.list | sort | uniq > fbox.uniq.list
perl /home/jfchen/FFproject/tools/bin/getidgff.pl -l fbox.uniq.list -g Glean.v1.4.chr.gff -o fbox.uniq.gff
perl /home/jfchen/FFproject/tools/findOverlap/findOverlap.pl fbox.uniq.mRNA.chr.table Gramene.mRNA.chr.table > fbox2Gramene

perl strand.pl fbox.uniq.gff ../Gramenev1.4/Gramene.chr.gff fbox2Gramene > fbox2Gramene.overlap
perl strand.pl Glean.v1.4.chr.gff ../Gramenev1.4/Gramene.chr.gff Glean2Gramene > Glean2Gramene.overlap


