perl single_family.pl --go ../input/Sb.final.pep.iprscan.gene.GO --species Sb > single_family_Sb.txt
perl single_family.pl --go ../input/Gramene.final.v1.4.pep.iprscan.gene.GO --species Ob > single_family_Ob.txt
perl single_family.pl --go ../input/tigr.all.final.iprscan.gene.GO --species Os > single_family_OS.txt

perl single_family2.pl --go ../input/Sb.final.pep.iprscan.gene.ipr --species Sb > single_family_Sb2.txt
perl single_family2.pl --go ../input/Gramene.final.v1.4.pep.iprscan.gene.ipr --species Ob > single_family_Ob2.txt
perl single_family2.pl --go ../input/tigr.all.final.iprscan.gene.ipr --species Os > single_family_OS2.txt

