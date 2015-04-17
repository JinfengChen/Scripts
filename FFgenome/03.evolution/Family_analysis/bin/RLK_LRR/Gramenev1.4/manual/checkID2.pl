#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"id1:s","id2:s","help");


my $help=<<USAGE;
Check if the gene id get by pfam is consistant with other publication
perl $0 --id1 --id2
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $refpfam=readtable($opt{id1});
my $refpub =readtable($opt{id2});
my ($miss,$add);
open OUT, ">miss.list" or die "$!";
foreach (sort keys %$refpub){
    #print "$_\n";
    if (exists $refpfam->{$_}){
       #print "$_\n";
    }else{
       $miss++;
       print OUT "$_\t$refpub->{$_}\n";
    }
}
#print "---------------------\n";
close OUT;
open OUT, ">add.list" or die "$!";
foreach (sort keys %$refpfam){
    if (exists $refpub->{$_}){
       #print "$_\n";
    }else{
       $add++;
       print OUT "$_\t$refpfam->{$_}\n";
    }
}
close OUT;
print "Miss Number: $miss\nAdd Number: $add\n";


#############################################
sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    #$unit[0]=$1 if ($unit[0]=~/(.*)\.\d+/);
    $hash{$unit[0]}= $unit[1] =~/\w+/ ? $unit[1] : 1;
}
close IN;
return \%hash;
}

