#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"fasta:s","gff:s","help");


my $help=<<USAGE;
This program will help you to get all repeat element in gff format.
perl getRepeat.pl -fasta /share/raid12/chenjinfeng/FFgenome/genome/scaffold/Scaffold.seq -gff allTE.gff -element > ../input/repeat.fa &
USAGE


if ($opt{help}){
print "$help\n";
exit;
}

$/=">";
my %seq;
open IN, "$opt{fasta}" or die "$!";
while (<IN>){
    chomp $_;
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $head=shift @unit;
    my @array=split(" ",$head);
    $head=$array[0];
    my $sequence=join("",@unit);
    $seq{$head}=$sequence;
}
close IN;
$/="\n";

open IN, "$opt{gff}" or die "$!";;
while (<IN>){
    chomp $_;
    next if ($_=~/^#/ or $_=~/^$/);
    my @unit=split("\t",$_);
    if ($unit[2] =~ /$opt{element}/){
        $unit[8]=~/ID=(.*);Target=(.*);Class=(.*?);/;
        my $id=$1."#".$3;
        my $start=$unit[3]-1;
        my $len  =abs($unit[4]-$unit[3])+1;
        if (exists $seq{$unit[0]}){
           my $element=substr($seq{$unit[0]},$start,$len);
           print ">$id\n$element\n"; 
        }else{
           print "Sequence $unit[0] not found!!!\n";
        }
    }
}
close IN;






