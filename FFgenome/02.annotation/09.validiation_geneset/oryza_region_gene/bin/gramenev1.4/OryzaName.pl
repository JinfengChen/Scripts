#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"fasta:s","region:s","help");


my $help=<<USAGE;
Format fasta sequence name for Oryza.fa
perl $0 --fasta ../input/Oryza.fa --region Moc1

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

& getfastaseq($opt{fasta});

sub getfastaseq
{
$/=">";
my %hash;
my ($file)=@_;
my $pre=$1 if ($file=~/(.*)\.fa/ or $file=~/(.*)\.txt/);
my $count=0;
open IN,"$file" or die "$!";
open OUT, ">$pre.final.fa" or die "$!";
while (<IN>){
    next if (length $_ < 2 or $_ =~/^$/);
    $count++;
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my $name=$opt{region}.".".$count;
    my $anno=$1 if ($temp=~/\[protein=(.*?)\]/);
    my $head="$name $anno";
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    print OUT ">$head\n$seq\n";
    $hash{$head}=$seq;
}
close IN;
close OUT;
$/="\n";
return \%hash;
}
 
