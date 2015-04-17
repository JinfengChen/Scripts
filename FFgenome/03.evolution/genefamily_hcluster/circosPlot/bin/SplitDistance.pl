#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"syntid:s","distance:s","help");


my $help=<<USAGE;
split ortholog distance into synteny or non-synteny distance
perl $0 --syntid --distance
--syntid: gene id file preduced by synteny_pipeline
--distance:
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 


my $refsyn=readtable($opt{syntid});
& ort($refsyn,$opt{distance});
###################
sub ort
{
my ($ref,$distance)=@_;
my ($n1,$n2);
my $name=$1 if ($distance=~/(.*)\.distance.txt/);
my $syn="$name.syn.distance.txt";
my $nonsyn="$name.nonsyn.distance.txt";
open OUT1, ">$syn" or die "$!";
open OUT2, ">$nonsyn" or die "$!";
open IN, "$distance" or die "$!";
<IN>;
while(<IN>){
   chomp $_;
   next if ($_=~/^$/);
   my @unit=split("\t",$_);
   if (exists $ref->{$unit[1]}){
      $n1++;
      print OUT1 "$_\n";
   }else{
      $n2++;
      print OUT2 "$_\n";
   }
}
close IN;
close OUT1;
close OUT2;
print "$name\tSynteny\tNonSynteny\n";
print "$name\t$n1\t$n2\n";
}

sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $unit[0] = $1 if ($unit[0]=~/(.*)\_\w+/);
    $hash{$unit[$unit[0]]}=1;
}
close IN;
return \%hash;
}


