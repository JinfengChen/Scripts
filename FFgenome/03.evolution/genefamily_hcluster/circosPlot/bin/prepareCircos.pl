#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"ortholog:s","gff1:s","gff2:s","chrlen1:s","chrlen2:s","help");


my $help=<<USAGE;
perl $0 --ortholog ../input/OBRACH_TIGR6.distance.txt --gff1 ../input/OBa.chr.gff --gff2 ../input/tigr.all.final.gff --chrlen1 ../input/OBa.chrlen --chrlen2 ../input/tigr.chrlen > log 2> log2 &

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}
my $reforth=ortholog($opt{ortholog});
my $refgff1=parseGFF($opt{gff1});
my $refgff2=parseGFF($opt{gff2});
my $refchr1=chrlen($opt{chrlen1});
my $refchr2=chrlen($opt{chrlen2});


my %chr;
my $count;
open OUT, ">karyotype.txt" or die "$!";
my $color;
foreach my $chr (sort {$a <=> $b} keys %$refchr1){
    $count++;
    $color++;
    print OUT "chr - chr$count Ob$chr 0 $refchr1->{$chr} chr$color\n";
    $chr{"Ob$chr"}="chr$count";
}
my $color;
foreach my $chr (sort {$a <=> $b} keys %$refchr2){
    $count++;
    $color++;
    print OUT "chr - chr$count Os$chr 0 $refchr2->{$chr} chr$color\n";
    $chr{"Os$chr"}="chr$count";
}
close OUT;


my $link;
open OUT, ">links.txt" or die "$!";
foreach my $id (keys %$reforth){
   if (exists $refgff1->{$reforth->{$id}->[0]}){
      $link++;
      my $gene1=$reforth->{$id}->[0];
      my $chr1=$refgff1->{$gene1}->[0];
      my $chr1a=$chr{"Ob$chr1"};
      my $end1=$refgff1->{$gene1}->[1]+10;
      #print OUT "link_$link $chr1a $refgff1->{$gene1}->[1] $refgff1->{$gene1}->[2] color=$chr1a\n";
      #print OUT "link_$link $chr1a $refgff1->{$gene1}->[1] $end1 color=$chr1a\n";
      
      my $gene2=$reforth->{$id}->[1];
      my $chr2=$refgff2->{$gene2}->[0];
      my $chr2a=$chr{"Os$chr2"}; 
      my $end2=$refgff2->{$gene2}->[1]+10;
      #print OUT "link_$link $chr2a $refgff2->{$gene2}->[1] $refgff2->{$gene2}->[2] color=$chr2a\n";
      print OUT "link_$link $chr1a $refgff1->{$gene1}->[1] $end1 color=$chr2a\n";
      print OUT "link_$link $chr2a $refgff2->{$gene2}->[1] $end2 color=$chr2a\n";
 
   }else{
      print "$reforth->{$id}->[0]\t$reforth->{$id}->[1]\n";
   }
}
close OUT;

####################################
sub parseGFF
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
        $unit[0]=~s/chr0//;
        $unit[0]=~s/chr//;
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        $hash{$id}=[$seq,$unit[3],$unit[4]];
    }
}
close IN;
return \%hash;
}

sub chrlen
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $unit[0]=~s/chr0//;
    $unit[0]=~s/chr//;
    $hash{$unit[0]}=$unit[1];
}
close IN;
return \%hash;
}

sub ortholog
{
my ($file)=@_;
my %hash;
my $count;
open IN, "$file" or die "$!";
#<IN>;
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    $count++;
    my @unit=split("\t",$_);
    $hash{$count}=[$unit[0],$unit[1]];
}
close IN;
return \%hash;
}

