#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"id4pfam:s","id4pub:s","help");


my $help=<<USAGE;
Check if the gene id get by pfam is consistant with other publication
perl $0 --id4pfam --id4pub
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $refpfam=readtable($opt{id4pfam});
my $refpub =readtable($opt{id4pub});
my ($miss,$add);
foreach (keys %$refpub){
    #print "$_\n";
    if (exists $refpfam->{$_}){
       #print "$_\n";
    }else{
       $miss++;
       print "Miss:$_\n";
    }
}
print "---------------------\n";
foreach (keys %$refpfam){
    if (exists $refpub->{$_}){
       #print "$_\n";
    }else{
       $add++;
       print "Add:$_\n";
    }
}

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
    $unit[0]=$1 if ($unit[0]=~/(.*)\.\d+/);
    $hash{$unit[0]}=1;
}
close IN;
return \%hash;
}

