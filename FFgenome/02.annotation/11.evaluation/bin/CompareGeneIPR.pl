#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"iprscan:s","dir:s","project:s","help");


my $help=<<USAGE;
perl $0 --iprscan --dir --project
--dir: dir of type list
count the number of gene that annotated by iprscan in each type of gene
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $refiprscan=readiprscan($opt{iprscan});
my @file=glob("$opt{dir}/$opt{project}.*");

my $gene;
foreach my $file (@file){
print "$file\n";
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    my $id=$1 if ($unit[0]=~/(.*)\_\w+/);
    if (exists $refiprscan->{$id}){
       $hash{1}++;
    }else{
       $hash{2}++;
    }
    $hash{3}++;
    $gene++;
}
close IN;
print "Gene in this file:\t$hash{3}\n";
print "Annotated gene number:\t$hash{1}\n";
print "Unannotated gene number:\t$hash{2}\n";
}
print "Total Gene:\t$gene\n";

sub readiprscan{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    unless (exists $hash{$unit[0]}){
       $hash{$unit[0]}=[$unit[11],$unit[12]];
    }  
}
close IN;
return \%hash;
} 
