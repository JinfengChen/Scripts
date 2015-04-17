#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"rice:s","distance:s","gff:s","help");


my $help=<<USAGE;
Get related gene of listed gene in rice from distance table.
These gene will be add to final8_3_best, so we need to check if these gene already exist in final8_3_best.
perl $0 --rice --distance --gff
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $ref=readtable($opt{rice});
my $gff=parseGFFid($opt{gff});

my %temp;
open IN, "$opt{distance}" or die "$!";
while(<IN>){
   chomp $_;
   my @unit=split("\t",$_);
   if (exists $ref->{$unit[1]}){
      unless (exists $gff->{$unit[0]}){
        unless (exists $temp{$unit[0]}){ 
          print "$unit[0]\n";
          $temp{$unit[0]}=1;
        }
      }
   }
}
close IN;
#OBR_GLEAN_10030480_OBRACH
sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=1;
}
close IN;
return \%hash;
}
 
sub parseGFFid
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;   ##Scaffold
my $id;    ##ID for element
my $record;##all line for this record
my $index; ##key, Seq_id
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*?);/){
            $id=$1;
        }
        $hash{$id}=$1;
    }

}
close IN;
return \%hash;
}

