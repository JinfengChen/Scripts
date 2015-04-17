#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"list:s","list2:s","help");


my $help=<<USAGE;
check these gene that exists in list but absent in list2.
perl $0 --list --list2
--list,list2: block.table
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $ref1=readtable1($opt{list});
my $ref2=readtable1($opt{list2});
foreach (keys %$ref1){
   unless (exists $ref2->{$_}){
      print "$_\n";
   }
}

#OBR_GLEAN_10030480_OBRACH
sub readtable1
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my @gene=split(";",$unit[3]);
    foreach my $g (@gene){
       $hash{$g}=1;
    }
}
close IN;
return \%hash;
}


