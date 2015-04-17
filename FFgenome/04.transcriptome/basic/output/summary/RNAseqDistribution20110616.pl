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


my ($total,$junction)=bed2table($opt{bed},"$opt{project}.reads");
my $gff2bed="GFF2BED.pl";

`perl $gff2bed -gff $opt{gene} -feature mRNA -bed $opt{project}.mRNA.bed`;
`perl $gff2bed -gff $opt{gene} -feature CDS -bed $opt{project}.CDS.bed`;

`/home/jfchen/FFproject/tools/BEDTools/bin/intersectBed -u -a $opt{bed} -b $opt{project}.mRNA.bed > $opt{project}.mRNA.overlap`;
`/home/jfchen/FFproject/tools/BEDTools/bin/intersectBed -u -a $opt{bed} -b $opt{project}.CDS.bed > $opt{project}.exon.overlap`;
`/home/jfchen/FFproject/tools/BEDTools/bin/intersectBed -u -a $opt{bed} -b $opt{repeat} > $opt{project}.repeat.overlap`;

my $repeatread=`wc -l $opt{project}.repeat.overlap`;
my $mRNAread=`wc -l $opt{project}.mRNA.overlap`;
my $exonread=`wc -l $opt{project}.exon.overlap`;

chomp ($repeatread,$mRNAread,$exonread);
my $intronread=$mRNAread-$exonread;
my $intergenic=$total-$mRNAread-$repeatread;
my $exonreadnojun=$exonread-$junction;
open OUT, ">$opt{project}.summary" or die "$!";

#print OUT "Total reads: $total\nJunction reads: $junction\n";
#print OUT "Repeat reads: $repeatread\nmRNA reads: $mRNAread\nNonJunction Exon reads: $exonreadnojun\n";
#print OUT "Intron reads: $intronread\nOther Intergenic reads: $intergenic\n";

print OUT "Total reads: $total\n";
print OUT "mRNA reads: $mRNAread\n";
print OUT "Repeat reads: $repeatread\n";
print OUT "OtherIntergenic reads: $intergenic\n";
print OUT "Junction reads: $junction\n";
print OUT "NonJunctionExon reads: $exonreadnojun\n";
print OUT "Intron reads: $intronread\n";


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

