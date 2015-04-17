#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"dir:s","help");


my $help=<<USAGE;
Clean these gene that present in other type in unknown.
perl $0 --dir 

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my @file=glob("$opt{dir}/*.fa");
my %hash;
foreach my $file (@file){
   next if ($file=~/$opt{dir}\/Unknown\.fa$/);
   getfastaseq($file,\%hash);
}

cleanfastaseq("$opt{dir}/Unknown.fa",\%hash);


############################################3
sub getfastaseq
{
$/=">";
my ($file,$ref)=@_;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    #print "$head\n";
    $ref->{$head}=$seq;
}
close IN;
$/="\n";
}

sub cleanfastaseq
{
$/=">";
my ($file,$ref)=@_;
open IN,"$file" or die "$!";
open OUT, ">Unknown.clean.fa" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    #print "$head\n";
    unless (exists $ref->{$head}){
       print OUT ">$head\n$seq\n";
    }
}
close IN;
close OUT;
$/="\n";
}
 
