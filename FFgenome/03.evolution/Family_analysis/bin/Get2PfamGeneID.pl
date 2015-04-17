#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"iprscan:s","domain1:s","domain2:s","help");


my $help=<<USAGE;
Get get id of given two Pfam domain id
perl $0 --iprscan ../input/iprscan/OB.iprscan --domain1 PF00560 --domain2 PF00069

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

& getid ($opt{iprscan},$opt{domain1},$opt{domain2});

sub getid
{
my ($file,$domain1,$domain2)=@_;
my %hash1;
my %hash2;
my $species;
if ($file=~/input\/iprscan\/(.*)\.iprscan/){
      $species=$1;
      #print "$species\n";
      open IN, "$file" or die "$!";
         while(<IN>){
            chomp $_;
            my @unit=split("\t",$_);
            next unless ($unit[3] eq "HMMPfam");
            $hash1{$unit[0]}=$unit[12] if ($unit[4] eq $domain1); 
            $hash2{$unit[0]}=$unit[12] if ($unit[4] eq $domain2);
         }
      close IN;
}
open OUT, ">$domain1.$domain2.$species.id" or die "$!";
foreach(keys %hash1){
   #print $_,"\n";
   print OUT "$_\n" if exists $hash2{$_};
   #print "$_\t$refid->{$_}\n";
}
close OUT;
} 
