#!/usr/bin/perl
use Getopt::Long;
use strict;
use warnings;

my %opt;

GetOptions (\%opt,"gene:s","repeat:s","unalign:s","project:s","help");


my $help=<<USAGE;
perl $0 --gene gene.gff --repeat repeat.gff --unalign unalign.table --project BAC

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $te=parseGFFrepeat($opt{repeat});
& gff2table($opt{gene},"CDS","$opt{project}.CDS");
& gff2table($opt{repeat},"Transposon","$opt{project}.repeat");
`perl findOverlap.pl $opt{unalign} $opt{project}.CDS > $opt{project}.CDS.overlap`;
`perl findOverlap.pl $opt{unalign} $opt{project}.repeat > $opt{project}.repeat.overlap`;
my ($total,$cds)=sumCDS("$opt{project}.CDS.overlap");
my $repeat=sumrepeat("$opt{project}.repeat.overlap",$te);

print "Total unalign :$total\n";
my $cdsrate=$cds/$total;
print "CDS: $cds\t$cdsrate\n";
foreach my $t (keys %$repeat){
   my $rate=$repeat->{$t}/$total;
   print "$t\t$repeat->{$t}\t$rate\n";
}
###################
sub parseGFFrepeat
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[8]=~/ID=(.*?);.*Class=(.*?);/){
       my $id=$1;
       my $class=$2;
       if ($class=~/LTR/){
          $class="LTR";
       }elsif($class=~/DNA/){
          $class="DNA";
       }else{
          $class="otherTE";
       }
       $hash{$id}=$class
    }
}
close IN;
return \%hash;
}

sub sumCDS
{
my ($file)=@_;
my $unalignlen;
my $incdslen;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split(" ",$_);
    $unalignlen+=$unit[1];
    next if ($unit[3] == 0);
    my @array=split(",",$unit[4]);
    $incdslen+=$array[2];
}
close IN;
return ($unalignlen,$incdslen);
}

sub sumrepeat
{
my ($file,$te)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split(" ",$_);
    next if ($unit[3] == 0);
    my %temp;
    for(my $i=4;$i<@unit;$i++){
       my @array=split(",",$unit[$i]);
       $temp{$array[0]}=$array[2];
    } 
    my @id=sort {$temp{$b} <=> $temp{$a}} keys %temp;
    my $id=$id[0];
    my $len=$temp{$id};
    my $type=$te->{$id};
    #print "$unit[0]\t$unit[2]\t$id\t$type\t$len\n";
    $hash{$type}+=$len;
}
close IN;
return \%hash;
}



sub gff2table
{
my ($gff,$feature,$table)=@_;
open IN, "$gff" or die "$!";
open OUT, ">$table.unsort" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    if ($unit[2] eq $feature and ($unit[8] =~/ID=(.*?);/ or $unit[8] =~/Parent=(.*?);/)){
        print OUT "$unit[0]\t$1\t$unit[3]\t$unit[4]\n";
    }
}
close IN;
close OUT;
system ("msort -k 1,n3 $table.unsort > $table");
system ("rm $table.unsort");
}
