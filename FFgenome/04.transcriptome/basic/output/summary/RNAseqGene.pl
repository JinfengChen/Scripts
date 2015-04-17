#!/usr/bin/perl
use Getopt::Long;
use strict;

my %opt;
GetOptions (\%opt,"bed:s","gene:s","project:s","help");


my $help=<<USAGE;
Summary coverage and hit of solexa read on each gene
-bed: bed format of read alignment
-gene: gff format of gene annotation
Run: perl $0 -bed testsum.bed -gene testsum.gene.gff -project testsum
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}
my $coverageBED="coverageBed";
my $gff2bed="GFF2BED.pl";

my $refgff=splitGFF($opt{gene});
foreach(sort keys %$refgff){
    my $scaf=$_;
    open TEMP, ">temp.gff" or die "$!";
        print TEMP "$refgff->{$scaf}";
    close TEMP;
    `grep "$scaf" $opt{bed} > temp.bed`;
    system ("perl $gff2bed -gff temp.gff -feature CDS -bed temp.exon.bed");   
    system ("$coverageBED -a temp.bed -b temp.exon.bed >> $opt{project}.exon.coverage");
}


#system ("perl $gff2bed -gff $opt{gene} -feature CDS -bed $opt{project}.exon.bed");
#system ("$coverageBED -a $opt{bed} -b $opt{project}.exon.bed > $opt{project}.exon.coverage");


my $totalgene=0;
my $covergene=0;
my $exonreads=0;
my $refexon=sumcoverage("$opt{project}.exon.coverage");
open OUT, ">$opt{project}.exon.coverage.table" or die "$!";
print OUT "Gene\tGene Length\tCover Length\tCover rate\tHit reads\n";
foreach(keys %$refexon){
    $totalgene++;
    my $temparray=$refexon->{$_};
    if ($temparray->[1] > 0){
       $covergene++;
    }
    $exonreads+=$temparray->[2];
    my $coverate=$temparray->[1]/$temparray->[0];
    print OUT "$_\t$temparray->[0]\t$temparray->[1]\t$coverate\t$temparray->[2]\n";
    
}
close OUT;
print "Total gene: $totalgene\tCover gene: $covergene\tExon reads: $exonreads\n";

#################################################3

sub splitGFF
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    next if ($_ eq "" or $_ =~/^#/);
    my @unit=split("\t",$_);
    if (exists $hash{$unit[0]}){
       $hash{$unit[0]}.=$_;
    }else{
       $hash{$unit[0]}=$_;
    }

}
close IN;
return \%hash;
}

sub sumcoverage
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
     chomp $_;
     next if ($_ eq "");
     my @unit=split("\t",$_);
     if (exists $hash{$unit[3]}){
         my $temp=$hash{$unit[3]};
         $temp->[0]+=$unit[8];
         $temp->[1]+=$unit[7];
         $temp->[2]+=$unit[6];
         $hash{$unit[3]}=$temp;
     }else{
         $hash{$unit[3]}=[$unit[8],$unit[7],$unit[6]];
         #print "$hash{$unit[3]}->[0]\n";
     }   
}
close IN;
return \%hash;
} 
