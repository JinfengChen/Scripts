#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"help");


my $help=<<USAGE;
Summary TF number for all file in ../input/sum
perl $0 > TF.family
USAGE


if ($opt{help}){
    print "$help\n";
    exit();
} 
my %anno;
my %ref;
my @tf=glob("../input/sum/*_tf_family");
foreach my $f (@tf){
    if ($f=~/sum\/(.*?)_tf_family/){
       my $species=$1;
       $ref{$species}=readtable($f,\%anno);
    }
}
print "Type";
foreach my $s (sort keys %ref){
    print "\t$s";
}
print "\tAnnotation\n";
foreach my $a (keys %anno){
    print "$a";
    foreach my $s (sort keys %ref){
        my $hash=$ref{$s};
        my $number;
        if (exists $hash->{$a}){
           $number=@{$hash->{$a}};
        }else{
           $number=0;
        }
        print "\t$number";
    }
    print "\t$anno{$a}\n";
}


sub readtable
{
my ($file,$anno)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $anno->{$unit[1]}=$unit[2] unless (exists $anno->{$unit[1]});
    push (@{$hash{$unit[1]}},$unit[0]);
}
close IN;
return (\%hash);
}

