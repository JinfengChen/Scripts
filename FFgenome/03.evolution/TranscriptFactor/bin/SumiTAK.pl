#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"help");


my $help=<<USAGE;


USAGE


if ($opt{help}){
    print "$help\n";
    exit();
} 

my ($ref,$reftype)=readtable($ARGV[0]);
foreach my $type (sort keys %$ref){
    my $n= @{$ref->{$type}};
    print "$type\t$n\t$reftype->{$type}\n";
}

sub readtable
{
my ($file)=@_;
my %hash;
my %anno;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $anno{$unit[1]}=$unit[2] unless (exists $anno{$unit[1]});
    push (@{$hash{$unit[1]}},$unit[0]);
}
close IN;
return (\%hash,\%anno);
}

