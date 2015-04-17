#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"iprscan:s","inf:s","help");


my $help=<<USAGE;
For these ortholog gene influenced by TE insertion, we pick out the iprscan result of them for WEGO drawing.
perl $0 --iprscan rice.iprscan --inf orth_OBa_rice_all.ortholog.inf
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 



#mC+|TE-
#TE-|mC+
#TE-|mC-
my ($ref1,$ref2)=readtable($opt{inf});
my $iprscan=readiprscan($opt{iprscan});

open OUT1, ">wego.silence1.txt" or die "$!";
open OUT2, ">wego.silence2.txt" or die "$!";
foreach(keys %$ref1){
   if (exists $iprscan->{$_}){
      my $line=join("\n",@{$iprscan->{$_}});
      print OUT1 "$line\n";
   }else{
      print "REF1\t$_\n";
   }
}
foreach(keys %$ref2){
   if (exists $iprscan->{$_}){
      my $line=join("\n",@{$iprscan->{$_}});
      print OUT2 "$line\n";
   }else{
      print "REF2\t$_\n";
   }
}

close OUT1;
close OUT2;

###########################
sub readtable
{
my ($file)=@_;
my %hash1;
my %hash2;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    if ($_=~/mC\+\|TE\-/){
       $hash1{$unit[6]}=$unit[0];
    }elsif($_=~/TE\-\|mC/){
       $hash2{$unit[6]}=$unit[0];
    }
}
close IN;
return (\%hash1,\%hash2);
}

sub readiprscan
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    push (@{$hash{$unit[0]}},$_);
}
close IN;
return \%hash;
}
