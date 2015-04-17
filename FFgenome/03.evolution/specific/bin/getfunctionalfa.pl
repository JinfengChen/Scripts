#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"id:s","seq:s","project:s","help");


my $help=<<USAGE;
get funtional fasta from all protein or cds.
--id: specific id
--seq: all fasta sequence
perl $0 -id OB.specific.id.txt -seq OB.fasta -project OB 
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $refid=getid($opt{id});
my $refseq=getfastaseq($opt{seq});

open OUT, ">$opt{project}.funtional.fa" or die "$!";
foreach(keys %$refseq){
   unless (exists $refid->{$_}){
      print OUT ">$_\n$refseq->{$_}\n";
   }
}


#########################
sub getid
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!"; 
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    unless (exists $hash{$_}){
        $hash{$_}=1;
    }
}
close IN;
return \%hash;
}

sub getfastaseq
{
$/=">";
my %hash;
my ($file)=@_;
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
    $hash{$head}=$seq;
}
$/="\n";
return \%hash;
}

