#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"tandem:s","align:s","help");


my $help=<<USAGE;
The script is designed to find tandem repeated gene in mcscan outfile (os_ob.aligns).

Run: perl findTandemInAlign.pl -tandem OB.tandem_repeat_gene.txt -align ob_os.aligns

USAGE

if ($opt{help} or keys %opt < 1){
   print "$help\n";
   exit();
}
my $refalign=parseAlign($opt{align});

open IN, "$opt{tandem}" or die "$!";
while(<IN>){
    chomp $_;
    my $line=$_;
    my @unit=split(" ",$_);
    foreach(@unit){
       #print "$_\n";
       if (exists $refalign->{$_}){
          print "$line\t$refalign->{$_}\n";
       }
    }
}
close IN;

sub parseAlign
{
my ($align)=@_;
my $flag;
my $num;
my %hash;
open IN, "$align" or die "$!";
while(<IN>){
   if ($_=~/^## Alignment (\d+):/){
      $num=$1;
      $flag=1;
   }elsif($flag and $_!=~/^#/){
      my @unit=split(" ",$_);
      #print "$unit[2]\t$unit[3]\t$unit[4]\n";
      $hash{$unit[2]}=$num;
      $hash{$unit[3]}=$num;
   } 
}
close IN;
return \%hash;
}
