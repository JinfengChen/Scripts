#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"hmmer3:s","domain:s","help");


my $help=<<USAGE;
Get gene id of given Pfam domain and record other Pfam domain if have
perl $0 --hmmer3 --domain

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

& getid ($opt{hmmer3},$opt{domain});

sub getid
{
my ($file,$domain)=@_;
my %hash;
my %anno;
my $species;
if ($file=~/\/(.*)\.hmmer3/ or $file=~/(.*)\.hmmer3/){
      $species=$1;
      open IN, "$file" or die "$!";
         while(<IN>){
            chomp $_;
            next if ($_=~/^#/ or $_=~/^$/);
            my @unit=split(" ",$_);
            next if ($unit[4] > 1);
            $unit[3]=$1 if ($unit[3]=~/(\w+?)\..*/);
            $anno{$unit[3]}=$unit[2] unless exists $anno{$unit[3]};
            $hash{$unit[0]}->{$unit[2]}=1;
         }
      close IN;
}
print "$domain.$species.id\n";
open OUT, ">$domain.$species.id" or die "$!";
my %type;
foreach(keys %hash){
   #print "$_\n";
   my $refhash=$hash{$_};
   my @pfam=sort {$a cmp $b} keys %$refhash;
   my $pfam=join(":",@pfam);
   #$type{$pfam}++ if ($pfam=~/$anno{$domain}/);
   push (@{$type{$pfam}},$_) if ($pfam=~/$anno{$domain}/);
   print OUT "$_\t$pfam\n" if ($pfam=~/$anno{$domain}/);
}
close OUT;
foreach(keys %type){
   my $num=@{$type{$_}};
   my $gene=join(",",@{$type{$_}});
   #print "$_\t$type{$_}\n";
   print "$_\t$num\t$gene\n";
}
} 
