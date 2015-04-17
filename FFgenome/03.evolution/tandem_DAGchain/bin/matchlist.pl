#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"gff:s","blasttable:s","project:s","help");


my $help=<<USAGE;
Create matchlist file for DAGCHAINER from gff and blasttable.

Run: perl $0 -gff --blasttable
-gff: GFF file 
-blasttable
-project:
USAGE

if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 


my $refgff=gff($opt{gff});

open OUT, ">$opt{project}.matchlist" or die "$!";
open IN, "$opt{blasttable}" or die "$!";
while(<IN>){
   chomp $_;
   next if ($_ =~ /^$/);
   my @unit=split("\t",$_);
   next unless (exists $refgff->{$unit[0]});
   print OUT "$refgff->{$unit[0]}\t$refgff->{$unit[4]}\t$unit[13]\n";
}
close IN;
close OUT;

sub gff
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    if ($unit[2] eq "mRNA"){
         my $start=$unit[3]; 
         my $end  =$unit[4];
         $unit[0]=~s/chr0//g;
         $unit[0]=~s/chr//g;
         my $chr =$unit[0];
         my $gene;
         if ($unit[8] =~/ID=(.*);/){
             $gene=$1;
         }
         $hash{$gene}="$chr\t$gene\t$start\t$end";
    }
}
close IN;
return \%hash;
}
