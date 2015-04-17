#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"geneid:s","ortholog:s","solor:s","hcluster:s","gff:s","project:s","help");


my $help=<<USAGE;
Find ortholog of interested gene id in list.
perl $0 --geneid ../input/flowering.id --ortholog ../input/OBRACH_TIGR6.distance.txt --solor ../input/all_vs_all.blast.m8.solar --hcluster ../input/all_vs_all.blast.m8.solar.forHC.hcluster --gff ../input/tigr.all.final.gff --project FlowerGene > log 2> log2 &
 
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

our $suffix="TIGR6";
my $refgene=readgene($opt{geneid},"0");
my $refdist=readdistance($opt{ortholog},"1");
my $reffamily=family($opt{hcluster});
my $refstart=genestart($opt{gff});

my ($total,$find,$no);
open OUT, ">$opt{project}.ortholog" or die "$!";
print OUT "Gene\tName\tPosition\tFamilyID\tOverlap\t$refdist->{Gene2}\n";
foreach $gene (keys %$refgene){
   $total++;
   my $family=$reffamily->{$gene};
   my $start=$refstart->{$gene};
   #my $gene1= $gene=~/(.*)\.\d+/ ? $1 : $gene;
   if (exists $refdist->{$gene}){
      $find++;
      my @unit=split("\t",$refdist->{$gene});
      #my $family=$reffamily->{$gene};
      #my $start=$refstart->{$gene};
      my $solor=`grep "$unit[0].*$unit[1]" $opt{solor}`;
      chomp $solor;
      next if ($solor eq "");
      my @array=split("\t",$solor);
      #print "$array[5]\t$array[6]\n";
      my $overlap=($array[8]-$array[7]+1)/$array[6]; ### overlap percentage in rice gene
      $overlap=int ($overlap*100);
      print OUT "$gene\t$refgene->{$gene}\t$start\t$family\t$overlap\t$refdist->{$gene}\n";
   }else{
      $no++;
      print "$gene\t$refgene->{$gene}\t$start\t$family\tNA\n";
   }
}
close OUT;

open OUT, ">$opt{project}.Summary" or die "$!";
print OUT "Total Gene:\t$total\n";
print OUT "Find Gene:\t$find\n";
print OUT "Miss Gene:\t$no\n";
close OUT;

##########################


sub readdistance
{
my ($file,$col)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    my $line=$_;
    my $gene= $unit[$col]=~/(.*)\.\d+/ ? $1 : $unit[$col];
    $hash{$gene}=$line;
}
close IN;
return \%hash;
}

sub readgene
{
my ($file,$col)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    #my $gene=$unit[$col];
    my $gene= $unit[$col]=~/(.*)\.\d+/ ? $1 : $unit[$col];
    $hash{$gene}=$unit[1];
}
close IN;
return \%hash;
}


sub family
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ =~/^$/);
    my @unit=split("\t",$_);
    my @gene=split(",",$unit[6]);
    foreach my $gene (@gene){
       if ($gene=~/(.*)\.\d+\_$suffix/){
          $hash{$1}=$unit[0];
       }
    }
}
close IN;
return \%hash;
}

sub genestart
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;   ##Scaffold
my $id;    ##ID for element
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*)\.\d+;/ or $unit[8] =~/ID=(.*)\.\d+/){
            $id=$1;
        }
        $hash{$id}=$unit[3]; 
    }
}
close IN;
return \%hash;
}


 
