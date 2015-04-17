#!/usr/bin/perl

use FindBin qw ($Bin);
use Getopt::Long;


GetOptions(\%opt,"refseq:s","start:s","end:s","repeat","gff:s","output:s","help:s");

my $help=<<USAGE;
Get gff for defined ref and coordinate
perl $0 -refseq chr01 --start 100 --end 1000000 --gff gene.gff --output sub.gff > log  
--repeat: if set to repeat, only filter by each line in gff file
USAGE

if (exists $opt{help} or keys %opt < 1){
   print $help;
   exit;
}

if ($opt{repeat}){
open IN, "$opt{gff}" or die "cat not open outfile: $opt{gff}\n";
open OUT, ">$opt{output}" or die "cat not open outfile: $opt{output}\n";
while(<IN>){
     chomp $_;
     my $line=$_;
     if ($line=~/^#/){
        print OUT "$line\n";
     }else{
        my @unit=split("\t",$line);
        if ($unit[0] eq $opt{refseq} and $unit[3] >= $opt{start} and $unit[4] <= $opt{end}){
           $unit[3] =$unit[3]-$opt{start}+1;
           $unit[4] =$unit[4]-$opt{start}+1;
           my $temp=join("\t",@unit);
           print OUT "$temp\n";
        }
     }
}
close IN;
close OUT;
}else{

my $refgff=parseGFF($opt{gff});
open OUT, ">$opt{output}" or die "cat not open outfile: $opt{output}\n";
 foreach my $g (keys %$refgff){
   my $gff=$refgff->{$g};
   my @line=split("\n",$gff);
   my @unit=split("\t",$line[0]);
   if ($unit[0] eq $opt{refseq} and $unit[3] >= $opt{start} and $unit[4] <= $opt{end}){
      foreach my $l (@line){
           my @array=split("\t",$l);
           $array[3] =$array[3]-$opt{start}+1;
           $array[4] =$array[4]-$opt{start}+1;
           my $temp=join("\t",@array);
           print OUT "$temp\n";

      }
   }
 }
close OUT;
}



sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;
my $id;    ##ID for element
my $record;##all line for this record
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
        $hash{$id}=$record;
    }elsif($unit[0] eq $seq and $unit[8] =~ /Parent=$id/){
        $hash{$id}.="$_\n";
    }

}
close IN;

return \%hash;
}


