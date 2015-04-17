#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"dir:s","help");


my $help=<<USAGE;
Classify gene into different type of homolog group (synort/nonsynort/nonort/unique).
perl $0 --dir 
--dir: directory of type files
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my @table=glob("$opt{dir}/*.table");

foreach my $file (@table){
   my $type=$1 if ($file=~/$opt{dir}\/\w+\.(.*)\.table$/);
   my $ref;
   if ($type=~/tandem/){
      $ref=readtable2($file);
   }else{
      $ref=readtable($file);
   }
   foreach my $g (keys %$ref){
      print "$g\t$type\n";
   }
}

########################
sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $unit[0]=$1 if ($unit[0]=~/(.*)\_\w+$/);
    $hash{$unit[0]}=1;
}
close IN;
return \%hash;
}

sub readtable2
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

