#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"fasta:s","id:s","output:s","help");


my $help=<<USAGE;
Assign MADS box gene id to TIGR id
perl $0 --fasta --id --output
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $refid=readtable($opt{id});
assignid($opt{fasta},$opt{output},$refid);


#############################################
sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
<IN>;
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    #print "$unit[0]\t$unit[1]\n";
    $hash{$unit[2]}="$unit[1]\_$unit[9]";
}
close IN;
return \%hash;
}

sub assignid
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
    my $locus=$1 if ($head=~/(.*)\.\d+/);
    #print "$locus\t$head\t$ref->{$locus}\n";
    if (exists $ref->{$head} or exists $ref->{$locus}){
        print OUT ">$head\_$ref->{$locus}\n$seq\n";
    }else{
        print OUT ">$head\n$seq\n";
    }
}
close IN;
close OUT;
$/="\n";
}

