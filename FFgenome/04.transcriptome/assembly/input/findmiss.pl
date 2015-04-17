#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"fasta:s","gene:s","genome:s","help");


my $help=<<USAGE;
perl $0 --fasta --gene --genome
find transcript that mapped to genome but not present in gene
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $gene=readpsl($opt{gene});
my $genome=readpsl($opt{genome});

foreach my $g (keys %$genome){
    unless (exists $gene->{$g}){
        print "$g\n";
    }
}


#####################
sub readpsl
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while (<IN>){
    my $line=$_;
    chomp $line;
    my @unit=split("\t",$_);
    $hash{$unit[9]}=1;
}
close IN;
return \%hash;
}
