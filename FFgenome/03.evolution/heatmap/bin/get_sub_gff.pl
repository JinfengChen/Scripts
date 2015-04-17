#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"list:s","gff:s","output:s","help");


my $help=<<USAGE;

Run: perl get_sub_gff.pl -list ricepara.gene -gff rice.gene.gff -output rice.para.gff
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $refgff=parseGFF($opt{gff});

open IN, "$opt{list}" or die "$!";
open OUT, ">$opt{output}" or die "$!";
while(<IN>){
   chomp $_;
   next if ($_ eq "");
   if (exists $refgff->{$_}){
      print OUT "$refgff->{$_}";
   }
}
close OUT;
close IN;
###########################################

sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;
my $id;    ##ID for element
my $record;##all line for this record
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        $record="$_\n";
        $hash{$id}=$record;
    }elsif($unit[0] eq $seq and $unit[8] =~ /Parent=$id/){
        $hash{$id}.="$_\n";
    }

}
close IN;

return \%hash;
}

