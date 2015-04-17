#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"pfam:s","tandem:s","help");


my $help=<<USAGE;
Find how many of Pfam gene is in tandem repeat gene set.

Run: perl findPfamInTandem.pl -pfam PF00646_OS.pfam -tandem OS.tandem_repeat_gene.txt

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $reftandem=tandem($opt{tandem});

my $total;
my $tandem;
open OUT, ">pfam.tandem" or die "$!";
open IN, "$opt{pfam}" or die "$!";
while(<IN>){
   chomp $_;
   next if ($_ eq "");
   $total++;
   if (exists $reftandem->{$_}){
       $tandem++;
       print OUT "$_\t$reftandem->{$_}\n"
   }
}
close IN;
close OUT;

print "Total Pfam: $total\n";
print "Tandem Pfam: $tandem\n";

#############################
sub tandem
{
my ($file)=@_;
my %hash;
my $counter;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    $counter++;
    my $line=$_;
    my @unit=split(" ",$_);
    foreach(@unit){
       $hash{$_}=$counter;
    }
}
close IN;
return \%hash;
}
