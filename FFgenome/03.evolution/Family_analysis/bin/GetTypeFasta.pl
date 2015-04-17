#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"type:s","fasta:s","help");


my $help=<<USAGE;
perl $0 --type --fasta
--type: type file of Fbox classification
--fasta: protein sequence
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my %hash;
my $project=$1 if ($opt{type}=~/(.*)\.type/);
open IN, "$opt{type}" or die "$!";
while(<IN>){
   chomp $_;
   next if ($_=~/^$/);
   my @unit=split("\t",$_);
   my @gene=split(",",$unit[2]);
   foreach my $gene (@gene){
      $hash{$gene}=$unit[0];
   }
}
close IN;

& typeseq($project,\%hash);
### read fasta file and output fasta if id is found in list file
sub typeseq
{
my ($project,$refhash)=@_;
$/=">";
open IN,"$opt{fasta}" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp,2);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    if (exists $refhash->{$head}){
      my $head1=$head."_$refhash->{$head}";
      open OUT ,">>$project.$refhash->{$head}.fa" or die "$!"; 
           print OUT ">$head1\n$seq";
      close OUT;
    }
}
close IN;
$/="\n";
}
