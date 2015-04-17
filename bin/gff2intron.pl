#!/usr/bin/perl
use Getopt::Long;
use strict;
use warnings;

my %opt;

GetOptions (\%opt,"gff:s","help");


my $help=<<USAGE;
perl $0 --gff
Generate intron and intergenic for gff
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}


my $refgff=parseGFF($opt{gff});
foreach my $chr (sort keys %$refgff){
   my $refchr=$refgff->{$chr};
   my $lastgend=1;
   foreach my $g (sort {$refchr->{$a}->[0] <=> $refchr->{$b}->[0]} keys %$refchr){
       #print "$chr\t$g\n";
       my @unit=split("\n",$refchr->{$g}->[1]);
       my $mrna=shift @unit;
       my @temp=split("\t",$mrna);
       my $temps=$lastgend;
       my $tempe=$temp[3]-1;
       $lastgend=$temp[4]+1;
       $temp[2]="intergenic";
       $temp[3]=$temps;
       $temp[4]=$tempe;
       my $inter=join("\t",@temp);
       print "$inter\n$mrna\n";
       my @cds;
       foreach my $cds (@unit){
          my @array=split("\t",$cds);
          push (@cds,[@array]);
       }
       @cds=sort {$a->[3] <=> $b->[3]} @cds;
       my $lasteend=$cds[0][4]+1;
       my $first=join("\t",@{$cds[0]});
       print "$first\n";
       for(my $i=1;$i<@cds;$i++){
           my $exon=join("\t",@{$cds[$i]});
           my @tempe=split("\t",$exon);
           my $tempis=$lasteend;
           my $tempie=$tempe[3]-1;
           $lasteend=$tempe[4]+1;
           $tempe[2]="intron";
           $tempe[3]=$tempis;
           $tempe[4]=$tempie;
           my $intron=join("\t",@tempe);
           print "$intron\n$exon\n";
       }
   }
}

=pod
my %hash;
open IN, "$opt{aligntable}" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    push @{$hash{$unit[0]}},[$unit[1],$unit[2],$unit[3]];
}
close IN;

foreach my $chr (sort keys %hash){
    my @array=@{$hash{$chr}};
    my $last=0;
    foreach my $block (sort {$a->[1] <=> $b->[1]} @array){
       my $end=$block->[1]-1;
       my $start=$last;
       $last=$block->[2]+1;
       print "$chr\t$block->[0]inter\t$start\t$end\n";
    }
}
=cut


sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;   ##Scaffold
my $id;    ##ID for element
my $record;##all line for this record
my $index; ##key, Seq_id
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*?);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        $record="$_\n";
        $hash{$seq}{$id}=[$unit[3],$record];
    }elsif($unit[2]=~/CDS/ and $unit[0] eq $seq and $unit[8] =~ /Parent=$id/){
        $hash{$seq}{$id}->[1].="$_\n";
    }

}
close IN;
return \%hash;
}

