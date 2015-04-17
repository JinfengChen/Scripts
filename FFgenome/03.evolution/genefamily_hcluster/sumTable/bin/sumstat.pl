#!/usr/bin/perl
use Getopt::Long;
use warnings;
use strict;

my %opt;

GetOptions (\%opt,"table:s","cluster:s","help");


my $help=<<USAGE;
perl $0 --table all_vs_all.blast.m8.solar.forHC.hcluster.stat --cluster all_vs_all.blast.m8.solar.forHC.hcluster
--table: stat table 
--cluster:
USAGE

$opt{table} ||= "all_vs_all.blast.m8.solar.forHC.hcluster.stat";
my %name=(
     "IRGSP" => "O. sativa",
     "TIGR6" => "O. sativa",
     "OBRACH" => "O. brachyantha",
     "BRDDI" => "B. distachyon",
     "BRADI" => "B. distachyon",
     "ZMAYS" => "Z. mays",
     "SBICO" => "S. bicolor",
     "POPRR"    => "P. trichocarpa"
);

my %total=(
     #"IRGSP" => "33265",
     "IRGSP" => "43230", # IRGSP with predicted gene
     "TIGR6" => "41342",
     #"OBRACH" => "30955",
     #"OBRACH" => "31059", #Gramene gene builder
     #"OBRACH" => "24915", #nr Gramene gene builder
     "OBRACH" => "32041", #final7 Gramene gene builder 34627
     "BRDDI" => "25532",
     "BRADI" => "25532",
     "ZMAYS" => "32540",
     "SBICO" => "34496",
     "POPRR"    => "45778"
);



if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

##column3: OBRACH
##column4: TIGR6
##column5: SBICO
my $cluster=readtable($opt{cluster});
my %summary;
open LOSE1, ">TIGR6_lose.gene" or die "$!";
open LOSE2, ">OBRACH_lose.gene" or die "$!";
open IN, "$opt{table}" or die "$!";
<IN>;
while(<IN>){
   chomp $_;
   next if ($_=~/^$/);
   my @unit=split("\t",$_);
   if ($unit[2] > 0 and $unit[3] ==0 and $unit[4] == 0){
      $summary{OBRACH_specific}[0]++;
      $summary{OBRACH_specific}[1]+=$unit[2];## number of specific gene
   }elsif($unit[2] == 0 and $unit[3] > 0 and $unit[4] == 0){
      $summary{TIGR6_specific}[0]++;
      $summary{TIGR6_specific}[1]+=$unit[3];
   }elsif($unit[2] == 0 and $unit[3] == 0 and $unit[4] > 0){
      $summary{SBICO_specific}[0]++;
      $summary{SBICO_specific}[1]+=$unit[4];
   }elsif($unit[2] > 0 and $unit[3] > 0 and $unit[4] == 0){
      $summary{ORYZA_specific}[0]++;
      $summary{ORYZA_specific}[1]+=$unit[2];
      $summary{ORYZA_specific}[2]+=$unit[3];
   }elsif($unit[2] > 0 and $unit[3] == 0 and $unit[4] > 0){
      writegene(\*LOSE1,"OBRACH",$cluster,$unit[0]);
      $summary{TIGR6_lose}[0]++;
      $summary{TIGR6_lose}[1]+=$unit[2]; ## number of gene another oryza species
   }elsif($unit[2] == 0 and $unit[3] > 0 and $unit[4] > 0){
      writegene(\*LOSE2,"TIGR6",$cluster,$unit[0]);
      $summary{OBRACH_lose}[0]++;
      $summary{OBRACH_lose}[1]+=$unit[3];
   }else{
      if ($unit[2] > $unit[3]){
         $summary{OBRACH_more}[0]++;
         $summary{OBRACH_more}[1]+=$unit[2];
         $summary{OBRACH_more}[2]+=$unit[3];
      }elsif($unit[2] < $unit[3]){
         $summary{TIGR6_more}[0]++;
         $summary{TIGR6_more}[1]+=$unit[2];
         $summary{TIGR6_more}[2]+=$unit[3];
      }else{
         $summary{OBRACH_TIGR6_SAME}[0]++;
         $summary{OBRACH_TIGR6_SAME}[1]+=$unit[2];## number of same gene
      }
   }
}
close IN;
close LOSE1;
close LOSE2;
#printhash(\%name);
#printhash(\%total);
printhash2(\%summary);


###################
sub writegene
{
my ($fh,$species,$cluster,$id)=@_;
foreach my $g (@{$cluster->{$id}}){
   print $fh "$g\t$id\n" if ($g=~/$species/);
}
}

sub printhash
{
my ($hash)=@_;
foreach (sort keys %$hash){
   print "$_\t$hash->{$_}\n";
}
}
##
sub printhash2
{
my ($hash)=@_;
foreach (sort keys %$hash){
   print "$_\t@{$hash->{$_}}\n";
}
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
    my @gene=split(",",$unit[6]);
    $hash{$unit[0]}=\@gene;
}
close IN;
return \%hash;
}

