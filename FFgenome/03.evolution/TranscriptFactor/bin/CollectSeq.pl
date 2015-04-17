#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"fasta:s","help");


my $help=<<USAGE;
Collect sequence for each type of TF
perl $0 --fasta
--fasta: pep file that merge all sequence
USAGE


if ($opt{help}){
    print "$help\n";
    exit();
} 
`mkdir ../input/sequence`; 
my $seq=getfastaseq($opt{fasta});
my @tf=glob("../input/sum/*_tf_family");
foreach my $f (@tf){
    if ($f=~/sum\/(.*?)_tf_family/){
       writeseq($f,$seq);
    }
}


sub writeseq
{
my ($file,$seq)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $tf=$unit[1];
    $tf=~s/-/_/g;
    $tf=~s/\//_/g;
    $tf=~s/ /_/g;
    open OUT, ">>../input/sequence/$tf.fasta" or die "$!";
         print OUT ">$unit[0]\n$seq->{$unit[0]}\n";
    close OUT;
}
close IN;
return (\%hash);
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
close IN;
$/="\n";
return \%hash;
}

