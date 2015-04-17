#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"hmmer3:s","domain1:s","domain2:s","help");


my $help=<<USAGE;
Get get id of given two Pfam domain id
perl $0 --hmmer3 tigr.hmmer3 --domain1 PF00560 --domain2 PF00069

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

& getid ($opt{hmmer3},$opt{domain1},$opt{domain2});

sub getid
{
my ($file,$domain1,$domain2)=@_;
my %hash1;
my %hash2;
my %anno;
my $species;
if ($file=~/\/(.*)\.hmmer3/ or $file=~/(.*)\.hmmer3/){
      $species=$1;
      #print "$species\n";
      open IN, "$file" or die "$!";
         while(<IN>){
            chomp $_;
            next if ($_=~/^#/ or $_=~/^$/);
            my @unit=split(" ",$_);
            next if ($unit[7] > 1);
            $unit[3]=$1 if ($unit[3]=~/(\w+?)\..*/);
            $anno{$unit[3]}=$unit[2] unless exists $anno{$unit[3]};
            
            $hash1{$unit[0]}=$unit[3] if ($unit[3] =~/$domain1/); 
            $hash2{$unit[0]}=$unit[3] if ($unit[3] =~/$domain2/);
         }
      close IN;
}
open OUT, ">$domain1.$domain2.$species.id" or die "$!";
foreach(sort keys %hash1){
   print $_,"\n";
   print OUT "$_\n" if exists $hash2{$_};
   #print "$_\t$refid->{$_}\n";
}
close OUT;
} 
