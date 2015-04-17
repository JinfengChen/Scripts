#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"iprscan:s","gff:s","help");


my $help=<<USAGE;
perl $0 --iprscan --gff
--gff: used to count total gene number
Get top 50 iprscan and count
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}


my $genen=$opt{gff} ? gffgenenumber($opt{gff}) : 30000 ;

my ($hash,$ipr)=readiprscan($opt{iprscan});

my @array;
foreach my $term (keys %$hash){
   my $num=keys %{$hash->{$term}};
   my $percent=$num/$genen;
   push @array, [$term,$num,$percent,$ipr->{$term}];
}

@temp =sort {$b->[1] <=> $a->[1]} @array;
for (my $i=0;$i<=$#temp;$i++){
    print "$temp[$i][0]\t$temp[$i][1]\t$temp[$i][2]\t$temp[$i][3]\n";
}

#################
sub readiprscan{
my ($file)=@_;
my %hash;
my %ipr;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    next if ($unit[11] eq "NULL");
    $ipr{$unit[11]}=$unit[12];
    $hash{$unit[11]}{$unit[0]}=1;
}
close IN;
return (\%hash,\%ipr);
}

sub gffgenenumber
{
my ($gff)=@_;
my $count;
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
       $count++;
    }

}
close IN;
return $count;
}
 
