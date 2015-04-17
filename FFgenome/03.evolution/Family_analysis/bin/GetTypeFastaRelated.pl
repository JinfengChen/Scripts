#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"type:s","hcluster:s","fasta:s","ob:s","sb:s","os:s","help");


my $help=<<USAGE;
Using the information of rice from type, we get related genes of OB/SB from hcluster to draw gene family tree.
perl $0 --type --fasta
--type: type file of Fbox classification
--hcluster:
--fasta: protein sequence
--ob: ob protein
--sb: sb protein
--os: os protein
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $refob=getfastaseq($opt{ob});
my $refsb=getfastaseq($opt{sb});
my $refos=getfastaseq($opt{os});
my %hash;
my %typegene;
my $project=$1 if ($opt{type}=~/(.*)\.type/);
open IN, "$opt{type}" or die "$!";
while(<IN>){
   chomp $_;
   next if ($_=~/^$/);
   my @unit=split("\t",$_);
   my @gene=split(",",$unit[2]);
   $typegene{$unit[0]}=\@gene;
   #my $tag=getrelated (\@gene,$unit[0],$opt{hcluster},$refob,$refsb,$refos,\%hash);
   #print $tag,"\n";
   #my %temp;
   foreach my $gene (@gene){
      $hash{$gene}=$unit[0];
      #$temp{$gene}=$unit[0];
   }
   #my $tag=getrelated (\@gene,$unit[0],$opt{hcluster},$refob,$refsb,$refos,\%temp);
}
close IN;

foreach my $type (keys %typegene){
   #print "$type\n";
   my $tag=getrelated ($typegene{$type},$type,$opt{hcluster},$refob,$refsb,$refos,\%hash);
}


#& typeseq($project,\%hash);
### read fasta file and output fasta if id is found in list file
sub typeseq
{
my ($project,$refhash)=@_;
$/=">";
open IN,"$opt{fasta}" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp,2);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    if (exists $refhash->{$head}){
      my $head1=$head."_$refhash->{$head}";
      open OUT ,">>$project.$refhash->{$head}.fa" or die "$!"; 
           print OUT ">$head1\n$seq";
      close OUT;
    }
}
close IN;
$/="\n";
}


sub getrelated
{
my ($gene,$type,$cluster,$ob,$sb,$os,$refhash)=@_;

my ($gene2clu,$clu2gene)=parsecluster($cluster);
my %fam;
### how many cluster 
foreach my $g (@$gene){
    #print "$g\n";
    $fam{$gene2clu->{$g}}=1;
}
### foreach cluster get ob and sb gene
open OUT1, ">OB.$type.cluster.fa" or die "$!";
open OUT2, ">SB.$type.cluster.fa" or die "$!";
open OUT3, ">OS.$type.cluster.fa" or die "$!";
foreach my $f (keys %fam){
    my $mem=$clu2gene->{$f};
    my @member=split(",",$mem);
    foreach my $m (@member){
       #print "$m\n";
       if ($m=~/(.*)\_OB/){
          my $id=$1;
          print OUT1 ">$id\n$ob->{$id}\n";
       }elsif($m=~/(.*)\_SB/){
          my $id=$1;
          print OUT2 ">$id\n$sb->{$id}\n";
       }elsif($m=~/(.*)\_TIGR/){
          my $id=$1;
          my $new;
          if (exists $refhash->{$id}){
             $new=$id."_".$refhash->{$id};
             print OUT3 ">$new\n$os->{$id}\n";
          }else{
             print OUT3 ">$id\n$os->{$id}\n";
          }
       }
    }
}
close OUT1;
close OUT2;
close OUT3;
return 1;
}


sub parsecluster
{
my ($file)=@_;
my %hash;
my %cluster;
open IN1, "$file" or die "$!";
while(<IN1>){
   chomp $_;
   my @unit=split("\t",$_);
   my @gene=split(",",$unit[6]);
   $cluster{$unit[0]}=$unit[6];
   foreach my $gene (@gene){
      if ($gene=~/(.*)\_(\w+)$/){
         $gene=$1;
         $hash{$gene}=$unit[0];
      }
   }
}
close IN1;
return (\%hash,\%cluster);
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
$/="\n";
return \%hash;
}

