#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"gff1:s","gff2:s","blastm8:s","inversion:s","help");


my $help=<<USAGE;
Get position and summary the inversion.
perl $0 --gff1 --gff2 --inversion
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}
my $blast=getblastm8($opt{blastm8});
my $gff1=parseGFFpos($opt{gff1});
my $gff2=parseGFFpos($opt{gff2});
open IN, "$opt{inversion}" or die "$!";
while(<IN>){
     chomp $_;
     my @unit=split("\t",$_);
     #print "$unit[1]\n";
     if ($unit[1]=~/\w+/){
         my $len=$gff1->{$unit[1]}->[2]-$gff1->{$unit[0]}->[1]+1;
         my $gene1=countgene($gff1,$unit[0],$unit[1],"refgene.list");
         my $genen1=@$gene1;
         my $gene2=countgene($gff2,$unit[2],$unit[3],"qrygene.list");
         my $genen2=@$gene2;
         my $pair =pairgene($blast,$gene1,$gene2,"pairgene.list");
         print "$gff1->{$unit[0]}->[0]\t$unit[0]\t$unit[1]\t$unit[2]\t$unit[3]\t$gff1->{$unit[0]}->[1]\t$gff1->{$unit[1]}->[2]\t$gff2->{$unit[2]}->[1]\t$gff2->{$unit[3]}->[2]\t$genen1\t$genen2\t$pair\t$len\n";      
     }else{
         open OUT, ">>singlepair.list" or die "$!";
              print OUT "$unit[0]\t$unit[2]\n";
         close OUT;
         my $len=$gff1->{$unit[0]}->[2]-$gff1->{$unit[0]}->[1]+1;
         print "$gff1->{$unit[0]}->[0]\t$unit[0]\tNA\t$unit[2]\tNA\t$gff1->{$unit[0]}->[1]\t$gff1->{$unit[0]}->[2]\t$gff2->{$unit[2]}->[1]\t$gff2->{$unit[2]}->[2]\t1\t1\t1\t$len\n";
     }
}
close IN;

######
sub pairgene
{
my ($blast,$gene1,$gene2,$file)=@_;
my %pair;
open OUT, ">>$file" or die "$!";
for(my $i=0;$i<@$gene1;$i++){
   for(my $j=0;$j<@$gene2;$j++){
       my $g1=$gene1->[$i];
       my $g2=$gene2->[$j];
       #print "$g1\t$g2\n";
       if (exists $blast->{"$g1&$g2"}){
          $pair{$g1}=1;  
          print OUT "$g1\t$g2\n";
       }
   }
}
close OUT;
my $p=keys %pair;
return $p;
}

sub countgene
{
my ($gff,$gene1,$gene2,$file)=@_;
my $chr=$gff->{$gene1}->[0];
my $start=$gff->{$gene1}->[1] > $gff->{$gene2}->[1] ? $gff->{$gene2}->[1] : $gff->{$gene1}->[1];
my $end  =$gff->{$gene1}->[1] > $gff->{$gene2}->[1] ? $gff->{$gene1}->[2] : $gff->{$gene2}->[2];
my @gene;
open OUT, ">>$file" or die "$!";
foreach my $g (keys %$gff){
   next if ($gff->{$g}->[0] ne $chr);
   if ($gff->{$g}->[1] >= $start and $gff->{$g}->[2] <= $end){
      push (@gene,$g);
      print OUT "$g\n";
   }
}
close OUT;
return \@gene;
}


sub getblastm8
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $unit[0]=$1 if ($unit[0]=~/LOC_(.*)_\w+$/ or $unit[0]=~/(.*)_\w+$/);
    $unit[1]=$1 if ($unit[1]=~/LOC_(.*)_\w+$/ or $unit[1]=~/(.*)_\w+$/);
    $hash{"$unit[0]&$unit[1]"}=1;
}
close IN;
return \%hash;
}


sub parseGFFpos
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $id;    ##ID for element
my $chr;
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        if ($unit[8]=~/ID=(.*?);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
            $id=$1 if ($id=~/LOC_(.*)/);
            $chr=$1 if ($unit[0]=~/chr(\d+)/);
        }
        $hash{$id}=[$chr,$unit[3],$unit[4]];       
    }
}
close IN;
return \%hash;
}

 
