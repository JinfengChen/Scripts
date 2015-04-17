#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"iprscan:s","domain:s","help");


my $help=<<USAGE;
Get get id of given Pfam domain id
perl $0 --iprscan --domain

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

& getid ($opt{iprscan},$opt{domain});

sub getid
{
my ($file,$domain)=@_;
my %hash;
my $species;
if ($file=~/input\/iprscan\/(.*)\.iprscan/){
      $species=$1;
      #print "$species\n";
      open IN, "$file" or die "$!";
         while(<IN>){
            chomp $_;
            my @unit=split("\t",$_);
            next unless ($unit[3] eq "HMMPfam");
            next unless ($unit[4] eq $domain);
            $hash{$unit[0]}=$unit[12];   
         }
      close IN;
}
open OUT, ">$domain.$species.id" or die "$!";
foreach(keys %hash){
   print OUT "$_\n";
   #print "$_\t$refid->{$_}\n";
}
close OUT;
} 
