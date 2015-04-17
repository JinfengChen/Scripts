#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"blastm8:s","query:s","target:s","help");


my $help=<<USAGE;
convert blastm8 to blasttable
perl $0 --query --target --blastm8
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $qrylen=getfastalen($opt{query});
my $reflen=getfastalen($opt{target});
my $title=$1 if ($opt{blastm8}=~/(.*)\.blast\.m8/);
convert($opt{blastm8},$qrylen,$reflen,$title);



##################
sub convert
{
my ($file,$qrylen,$reflen,$title)=@_;
open OUT, ">$title.blasttable" or die "$!";
open IN, "$file" or die "$!";
while(<IN>){
   chomp $_;
   my @unit=split("\t",$_);
   print OUT "$unit[0]\t$qrylen->{$unit[0]}\t$unit[6]\t$unit[7]\t$unit[1]\t$reflen->{$unit[1]}\t$unit[8]\t$unit[9]\t$unit[2]\t--\t$unit[5]\t$unit[3]\t$unit[11]\t$unit[10]\t--\t--\n";   
}
close IN;
close OUT;
}


sub getfastalen
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
    $hash{$head}= length $seq;
}
close IN;
$/="\n";
return \%hash;
}

