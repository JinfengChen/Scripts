#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"blast:s","position:s","help");


my $help=<<USAGE;
Get closest homolog of gene in spe1 in other two species by compare the score of blast.
perl $0 --blast --position
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}
my $pos=readtable($opt{position});
my ($id,$homo1,$homo2)=blastpair($opt{blast});
my $count;
foreach my $g ( keys %$id){
   $count++;
   my $pos0= $pos->{$g};
   my $hit1= $homo1->{$g}->[0] ? $homo1->{$g}->[0] : "NA";
   my $pos1= $pos->{$hit1} ? $pos->{$hit1} : "NA";
   my $hit2= $homo2->{$g}->[0] ? $homo2->{$g}->[0] : "NA";
   my $pos2= $pos->{$hit2} ? $pos->{$hit2} : "NA";
   print "$g\t$pos0\t$hit1\t$pos1\t$hit2\t$pos2\t$count\tNA\n";
}
#########
sub blastpair
{
my ($blast)=@_;
my $spe1="Os";
my $spe2="Sb";
my %id;
my %homolog1; ## for rice
my %homolog2; ## for sorghum
open IN, "$blast" or die "$!";
<IN>;
while(<IN>){
   my @unit=split("\t",$_);
   if ($unit[13] > 1e-10){
       next;
   }
   $id{$unit[0]}=1;
   if ($unit[4]=~/$spe1/){
          $homolog1{$unit[0]} = $unit[12] > $homolog1{$unit[0]}->[1] ? [$unit[4],$unit[12]] : $homolog1{$unit[0]};
   }
   if ($unit[4]=~/$spe2/){
          $homolog2{$unit[0]} = $unit[12] > $homolog2{$unit[0]}->[1] ? [$unit[4],$unit[12]] : $homolog2{$unit[0]};
   }
}
close IN;
return (\%id,\%homolog1,\%homolog2);
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
    $hash{$unit[0]}=$unit[2];
}
close IN;
return \%hash;
}
 
