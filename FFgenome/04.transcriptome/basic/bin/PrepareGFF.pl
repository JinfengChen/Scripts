#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"gff:s","help");


my $help=<<USAGE;


USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 
my $head;
if ($opt{gff}=~/(.*)\.gff/){
   $head=$1;
}
my $refgff=parseGFF($opt{gff});

my $exon;
open OUT, ">$head.tophat.gff" or die "$!";
foreach(keys %$refgff){
    my $id=$_;
    my @line=split("\n",$refgff->{$_});
    foreach(@line){
       my @unit=split("\t",$_);
       if ($unit[2] eq "mRNA"){
          
          $unit[8]="ID=$id"."_T01;Parent=$id;";
          my $line2=join("\t",@unit);
          $unit[2]="gene";
          $unit[8]="ID=$id;";
          my $line1=join("\t",@unit);
          print OUT "$line1\n$line2\n";
       }elsif($unit[2] eq "CDS"){
          $exon++;
          $unit[2]="exon";
          $unit[8]="ID=exon$exon;Parent=$id\_T01;";
          my $line1=join("\t",@unit);
          print OUT "$line1\n";
       }else{
          print OUT "$_\n";
       }
    }
}
close OUT;

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
        if ($unit[8]=~/ID=(.*);/ or $unit[8] =~/ID=(.*)/){
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


