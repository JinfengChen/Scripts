#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"bed:s","gene","repeat","methy:s","project:s","help");


my $help=<<USAGE;
Run methylevelDZ.pl and drawmethyDZ.pl for gene feature.
--gene: select to deal with gene
--repeat: select to deal with repeat
--bed:  dir of bed file
--methy: dir of methylation bed file for CG,CHG,CHH
--project: project name
perl $0 -bed FF_gene_bed -methy FF_methylation_bed -gene -project FF
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

`mkdir $opt{project}` unless (-f $opt{project});
if ($opt{gene}){
     `perl methylevelDZ.pl -bed $opt{bed} -feature gene -type CG -methy $opt{methy} -project $opt{project} > CG.log`;
     `perl methylevelDZ.pl -bed $opt{bed} -feature gene -type CHG -methy $opt{methy} -project $opt{project} > CHG.log`;
     `perl methylevelDZ.pl -bed $opt{bed} -feature gene -type CHH -methy $opt{methy} -project $opt{project} > CHH.log`;
     my $cgsum="$opt{project}.gene.CG.level.sum";
     my $chgsum="$opt{project}.gene.CHG.level.sum";
     my $chhsum="$opt{project}.gene.CHH.level.sum";
     my $title="$opt{project}.Gene";
     `perl drawmethyDZ.pl -CG $cgsum -CHG $chgsum -CHH $chhsum -project $title > draw.log`;
     `mv $opt{project}.* $opt{project}`;
     `mv *.log $opt{project}`;
} 

my @rep=("REPEAT","RT","DNA","COPIA","GYPSY","MITE","CACTA","MUDR");
#my @rep=("MITE");
if ($opt{repeat}){
  for(my $i=0;$i<@rep;$i++){
    `perl methylevelDZ.pl -bed $opt{bed} -feature $rep[$i] -type CG -methy $opt{methy} -project $opt{project} > CG.log`;
    `perl methylevelDZ.pl -bed $opt{bed} -feature $rep[$i] -type CHG -methy $opt{methy} -project $opt{project} > CHG.log`;
    `perl methylevelDZ.pl -bed $opt{bed} -feature $rep[$i] -type CHH -methy $opt{methy} -project $opt{project} > CHH.log`;
     my $cgsum="$opt{project}.$rep[$i].CG.level.sum";
     my $chgsum="$opt{project}.$rep[$i].CHG.level.sum";
     my $chhsum="$opt{project}.$rep[$i].CHH.level.sum";
     my $title="$opt{project}.$rep[$i]";
     `perl drawmethyDZ.pl -CG $cgsum -CHG $chgsum -CHH $chhsum -project $title > draw.log`;
     `mv $opt{project}.* $opt{project}`;
     `mv *.log $opt{project}`;
  }
}




