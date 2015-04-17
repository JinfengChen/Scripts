#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"distance:s","spe1:s","spe2:s","help");


my $help=<<USAGE;
Summary ortholog and distance for two species from ortholog_paralog_pipeline.
perl $0 --distance ./ --spe1 OB -spe2 OS > log 2> log2 &
--distance: directory have all.*
--spe1: suffix for species 1, Os09t0109500-01_IRGSP is IRGSP
--spe2: suffix for species 2
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $orth="$opt{distance}/all.ortho";
my $kaks="$opt{distance}/all.KaKs";
my $iden="$opt{distance}/all.identity";
my $fdtv="$opt{distance}/all.4dtv";

my $refiden=getidentity($iden);
print "identity done\n";
my $reffdtv=get4dtv($fdtv);
print "4dtv done\n";
my $refkaks=getkaks($kaks);
print "kaks done\n";

print "$orth\n";
my %gene1;
my %gene2;
my %hash;
my $num;
#open OUT, ">$opt{spe1}2$opt{spe2}.distance.txt" or die "$!";
#print OUT "Gene1\tGene2\tCDS_identity\tprotein_identity\t4dtv_corrected\t4dtv_raw\tKa\tKs\tKaKs\n";
open IN, "$orth" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $pair=$unit[0];
   
    if ($pair=~/$opt{spe1}/ and $pair=~/$opt{spe2}/) {
       my @array=split("&",$pair);
       #print "$array[0]\t$array[1]\n";
       if ($pair=~/(.*)\_$opt{spe1}\&(.*)\_$opt{spe2}/){
          $hash{$pair}=[$1,$2,"$refiden->{$pair}\t$reffdtv->{$pair}\t$refkaks->{$pair}"];
          #print OUT "$1\t$2\t$refiden->{$pair}\t$reffdtv->{$pair}\t$refkaks->{$pair}\n";
          $gene1{$1}++;
          $gene2{$2}++;
       }elsif($pair=~/(.*)\_$opt{spe2}\&(.*)\_$opt{spe1}/){
          $hash{$pair}=[$2,$1,"$refiden->{$pair}\t$reffdtv->{$pair}\t$refkaks->{$pair}"];
          #print OUT "$2\t$1\t$refiden->{$pair}\t$reffdtv->{$pair}\t$refkaks->{$pair}\n";
          $gene1{$2}++;
          $gene2{$1}++;
       } 
       #print "$pair\n"; 
       $num++;
    }    
}
close IN;

open OUT, ">$opt{spe1}2$opt{spe2}.distance.txt" or die "$!";
print OUT "Gene1\tGene2\tOrthologType\tCDS_identity\tprotein_identity\t4dtv_corrected\t4dtv_raw\tKa\tKs\tKaKs\n";
foreach(keys %hash){
    my $pair=$_;
    if ($gene1{$hash{$pair}->[0]} == 1 and $gene2{$hash{$pair}->[1]} == 1){
       print OUT "$hash{$pair}->[0]\t$hash{$pair}->[1]\tortholog_one2one\t$hash{$pair}->[2]\n";
    }elsif($gene1{$hash{$pair}->[0]} > 1 and $gene2{$hash{$pair}->[1]} > 1){
       print OUT "$hash{$pair}->[0]\t$hash{$pair}->[1]\tortholog_many2many\t$hash{$pair}->[2]\n";
    }else{
       print OUT "$hash{$pair}->[0]\t$hash{$pair}->[1]\tortholog_one2many\t$hash{$pair}->[2]\n";
    }

}
close OUT;

my $num1=keys %gene1;
my $num2=keys %gene2;
print "Total Ortholog Pair:\t$num\n";
print "$opt{spe1} Gene Number:\t$num1\n";
print "$opt{spe2} Gene Number:\t$num2\n";


sub getidentity
{
my ($file)=@_;
my %hash;
open IN ,"$file" or die "$!";
while(<IN>){
   chomp $_;
   next if ($_=~/^$/);
   my @unit=split("\t",$_);
   my $pair=$unit[0];
   ### $unit[1] is cds identity $unit[2] is protein identity
   if ($pair=~/$opt{spe1}/ and $pair=~/$opt{spe2}/) {
      $hash{$pair}="$unit[1]\t$unit[2]";
   }
}
close IN;
return \%hash;
}

sub get4dtv
{
my ($file)=@_;
my %hash;
open IN ,"$file" or die "$!";
while(<IN>){
   chomp $_;
   next if ($_=~/^$/);
   my @unit=split("\t",$_);
   my $pair=$unit[0];
   ### $unit[1] is corrected 4dtv $unit[2] is raw 4dtv
   if ($pair=~/$opt{spe1}/ and $pair=~/$opt{spe2}/) {
      $hash{$pair}="$unit[1]\t$unit[2]";
   }
}
close IN;
return \%hash;
}


sub getkaks
{
my ($file)=@_;
my %hash;
open IN ,"$file" or die "$!";
while(<IN>){
   chomp $_;
   next if ($_=~/^$/);
   my @unit=split("\t",$_);
   my $pair=$unit[0];
   ### $unit[1] is Ka $unit[2] is Ks $unit[3] is Ka/Ks
   if ($pair=~/$opt{spe1}/ and $pair=~/$opt{spe2}/) {
       $hash{$pair}="$unit[2]\t$unit[3]\t$unit[4]";
   }
}
close IN; 
return \%hash; 
} 


