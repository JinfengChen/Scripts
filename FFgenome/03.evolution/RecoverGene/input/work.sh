awk -F"\t" '{print $1,$4,$5,$3,$6,$2}' super-scaffold.domain > super-scaffold.domain.bed
grep "mRNA" super-scaffold.gff > super-scaffold.mRNA.gff
perl ../bin/GFF2BED.pl --gff super-scaffold.gff --feature mRNA --bed super-scaffold.mRNA.bed

