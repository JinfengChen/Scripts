#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"cluster:s","help");


my $help=<<USAGE;
Read in the result of hcluster and pep sequence for each speices.
The output files will be sepcific gene for each species, which means these genes not present in any cluster.

perl $0 -cluster all_vs_all.blast.m8.solar.forHC.hcluster
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}
my %hash;
my %gene;
open IN, "$opt{cluster}" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    my @gene=split(",",$unit[6]);
    foreach(@gene){
       if ($_=~/(.*)\_(\w+)$/){
           if (exists $hash{$2}){
               my $reftemp=$hash{$2};
               $reftemp->{$1}=1;
               $hash{$2}=$reftemp;
           }else{
               my %temp;
               $temp{$1}=1;      
               $hash{$2}=\%temp;
           }
       }
    }    
}
close IN;

 
foreach(keys %hash){
   my $species=$_;
   my $iprfile="../input/$_.ipr";
   my $pepfile="../input/$_.pep";
   next unless (-f $iprfile);
   my $refipr=ipr($iprfile);
   my $refpep=getfastaid($pepfile);
   my $reftemp=$hash{$species};
   my @gene=keys %$reftemp;
   foreach (@gene){
      if (exists $refpep->{$_}){
         delete $refpep->{$_};
      }
   }   
   open OUT, ">$species.specific" or die "$!";
   foreach(keys %$refpep){ 
      my $ann;
      if (exists $refipr->{$_}){
         $ann=$refipr->{$_};
      }else{
         $ann="NA";
      }
      print OUT "$_\t$ann\n";
   }
   close OUT;
}


sub getfastaid
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
    $hash{$head}=1;
}
$/="\n";
return \%hash;
}



sub ipr
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    unless (exists $hash{$unit[0]}){
       $hash{$unit[0]}=$unit[2];
    } 
}
close IN;
return \%hash;
}


