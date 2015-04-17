#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"iprscan:s","domain:s","help");


my $help=<<USAGE;
Get get id of given Pfam domain and record other Pfam domain if have
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
my %anno;
my $species;
if ($file=~/input\/iprscan\/(.*)\.iprscan/){
      $species=$1;
      open IN, "$file" or die "$!";
         while(<IN>){
            chomp $_;
            my @unit=split("\t",$_);
            next unless ($unit[3] eq "HMMPfam");
            $anno{$unit[4]}=$unit[5] unless exists $anno{$unit[4]};
            $hash{$unit[0]}->{$unit[5]}=1;
         }
      close IN;
}
open OUT, ">$domain.$species.id" or die "$!";
my %type;
foreach(keys %hash){
   my $refhash=$hash{$_};
   my @pfam=sort {$a cmp $b} keys %$refhash;
   my $pfam=join(":",@pfam);
   $type{$pfam}++ if ($pfam=~/$anno{$domain}/);
   print OUT "$_\t$pfam\n" if ($pfam=~/$anno{$domain}/);
}
close OUT;
foreach(keys %type){
   print "$_\t$type{$_}\n";
}
} 
