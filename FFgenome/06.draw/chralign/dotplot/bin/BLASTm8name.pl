#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"blastm8:s","exclude:s","target:s","project:s","help");

my $help=<<USAGE;
Write for quote-align, do not used in for mcscan.
Run: perl $0 -blastm8 all_vs_all.blast.m8 -project 3way
-blast: blasttable file
-exclude: exclude species lines
-target:
-project: prefix for result blast file

USAGE

if ($opt{help} or keys %opt < 1){
   print "$help\n";
   exit();
}


my $filter_blast="/home/jfchen/FFproject/tools/mcscan0.8/filter_blast.py";

my $unfiltered=$opt{blastm8}.".formatname4$opt{project}";
& name($opt{blastm8},$unfiltered);


sub name
{
my ($file,$out)=@_;
open IN, "$file" or die "$!";
open OUT, ">$out" or die "$!";
while(<IN>){
   chomp $_;
   next if ($_=~/^$/);
   my @unit=split("\t",$_);
   next if ($unit[0]=~/$opt{exclude}/ or $unit[1]=~/$opt{exclude}/);
   next unless ($unit[0]=~/TIGR6/ and $unit[1]=~/$opt{target}/);
   $unit[0]=$1 if ($unit[0]=~/(.*)\_\w+$/);
   $unit[1]=$1 if ($unit[1]=~/(.*)\_\w+$/);
   my $line=join("\t",@unit);
   print OUT "$line\n";
}
close IN;
close OUT;
}


