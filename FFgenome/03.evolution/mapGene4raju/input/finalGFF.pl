#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"gff1:s","gff2:s","help");


my $help=<<USAGE;
perl $0 --gff1 final_v2.obra.gff --gff2 map.gff
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}
my $gffhead=GFFhead($opt{gff1});
my $gff=parseGFF($opt{gff2});


foreach my $g (keys %$gff){
   next if ($g=~/\_D\d+/);
   #print "$gff->{$g}";
   my @unit=split("\n",$gff->{$g});
   my $head=shift @unit;
   $gffhead->{$g}=~s/\r//;
   $gffhead->{$g}=~s/\n//;
   $head.="$gffhead->{$g}";
   print "$head\n";
   my $line=join("\n",@unit);
   print "$line\n";
}



###########################
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
    if ($unit[2]=~/mRNA/ or $unit[2]=~/predicted_gene/){
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*?);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        $record="$_\n";
        $index="$id";
        $hash{$index}=$record;
    }elsif($unit[2]=~/CDS/ and $unit[0] eq $seq and $unit[8] =~ /Parent=$id/){
        $hash{$index}.="$_\n";
    }

}
close IN;
return \%hash;
}

sub GFFhead
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/predicted_gene/){
        if ($unit[8]=~/ID=(.*?);(.*)/){
           $hash{$1}=$2; 
        }
    }

}
close IN;
return \%hash;
}
 
