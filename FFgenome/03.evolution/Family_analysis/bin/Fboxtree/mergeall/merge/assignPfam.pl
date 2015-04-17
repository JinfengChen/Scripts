#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"dir:s","help");


my $help=<<USAGE;
Assign pfam domain inf to gene
perl $0 --dir ./
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my @file=glob("$opt{dir}/*.fa");
my %hash;
getpfam("$opt{dir}/Pfam.fa",\%hash);

foreach my $file (@file){
   if ($file=~/$opt{dir}\/Pfam\.fa$/){
       next;
   }elsif($file=~/$opt{dir}\/(.*)\.fa$/){
       my $type=$1;
       my $out=$type.".pfam".".fa";
       assignpfam($file,$out,\%hash);
   }
}


############################################3
sub getpfam
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
    my $gene;
    my $type;
    if ($head=~/(.*)\_(\w+)$/){
        $gene=$1;
        $type=$head;
    }
    $ref->{$gene}=$type;
}
close IN;
$/="\n";
}
 
sub assignpfam
{
$/=">";
my ($file,$out,$ref)=@_;
open IN,"$file" or die "$!";
open OUT, ">$out" or die "$!";
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
    if (exists $ref->{$head}){
        print OUT ">$ref->{$head}\n$seq\n";
    }else{
        print OUT ">$head\n$seq\n";
    }
}
close IN;
close OUT;
$/="\n";
}

