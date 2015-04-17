#!/usr/bin/perl
use Getopt::Long;
use strict;

my %opt;
GetOptions (\%opt,"bed:s","gene:s","repeat:s","project:s","help");


my $help=<<USAGE;
The scripts summary the result of tophat and report how the RNAseq reads are distributed in each feature, exon/intron/TE/intergenic. 
Note: tophat should be run with -G, --no-novel-juncs, so any read that have a map length > 75 is junctions reads.


-bed: accepted_hit.bed, tophat result in bed format
-gene:   gene anntation file in gff format
-repeat: repeat anntation file in gff format
Run: perl $0 -bed accepted_hit.bed -gene FF.gene.gff -repeat FF.repeat.gff -project test
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $findoverlap ="findOverlap.pl";
my $gff2bed="GFF2BED.pl";


## do not use awk for bed alignment, write a sub function to count junction read and convert only intact read into table.
#`awk '{print \$1,\$4,\$2+1,\$3}' $opt{bed} > $opt{project}.reads`;
my ($total,$junction)=bed2table($opt{bed},"$opt{project}.reads");

& gff2table($opt{repeat},"Transposon","$opt{project}.repeat");
& gff2table($opt{gene},"mRNA","$opt{project}.mRNA");
& gff2table($opt{gene},"CDS","$opt{project}.exon");


#`awk '{if(\$3~/Transposon/){print \$1,\$2,\$4,\$5}}' $opt{repeat} > $opt{project}.repeat`;
#`awk '{if(\$3~/mRNA/){print \$1,\$2,\$4,\$5}}' $opt{gene} > $opt{project}.mRNA`;
#`awk '{if(\$3~/CDS/){print \$1,\$2,\$4,\$5}}' $opt{gene} > $opt{project}.exon`;

`perl $findoverlap $opt{project}.reads $opt{project}.repeat > $opt{project}.repeat.overlap`;
`perl $findoverlap $opt{project}.reads $opt{project}.mRNA > $opt{project}.mRNA.overlap`;
`perl $findoverlap $opt{project}.reads $opt{project}.exon > $opt{project}.exon.overlap`;

my $repeatread=`awk '\$4 > 0' $opt{project}.repeat.overlap | wc -l`;
my $mRNAread=`awk '\$4 > 0' $opt{project}.mRNA.overlap | wc -l`;
my $exonread=`awk '\$4 > 0' $opt{project}.exon.overlap | wc -l`;

chomp ($repeatread,$mRNAread,$exonread);
my $intronread=$mRNAread-$exonread;
my $intergenic=$total-$mRNAread-$junction-$repeatread;

open OUT, ">$opt{project}.summary" or die "$!";

print OUT "Total reads: $total\nJunction reads: $junction\n";
print OUT "Repeat reads: $repeatread\nmRNA reads: $mRNAread\nExon reads: $exonread\n";
print OUT "Intron reads: $intronread\nIntergenic reads: $intergenic\n";

close OUT;
#########################################################
sub bed2table
{
my ($bed,$read)=@_;
my $total=0;
my $junction=0;
open IN, "$bed" or die "$!";
open OUT, ">$read" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    $total++;
    my @unit=split("\t",$_);
    my $len=$unit[2]-$unit[1]; ### BED start position is 0-based, end is 1-based so the length=end-start
    if ($len > 75){
       $junction++;
    }else{
       my $start=$unit[1]+1;
       print OUT "$unit[0]\t$unit[3]\t$start\t$unit[2]\n";
    }
}
close IN; 
close OUT;
return ($total,$junction);
}

sub gff2table
{
my ($gff,$feature,$table)=@_;
open IN, "$gff" or die "$!";
open OUT, ">$table.unsort" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    if ($unit[2] eq $feature and ($unit[8] =~/ID=(.*)?;/ or $unit[8] =~/Parent=(.*)?;/)){
        print OUT "$unit[0]\t$1\t$unit[3]\t$unit[4]\n"; 
    }
}
close IN;
close OUT;
system ("msort -k 1,n3 $table.unsort > $table");
system ("rm $table.unsort");
}
