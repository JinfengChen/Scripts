#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"gff1:s","gff2:s","help");


my $help=<<USAGE;
perl $0 --gff1 --gff2

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}
my $gff1=parseGFF($opt{gff1});
my $gff2=parseGFF($opt{gff2});

my $exonnumber1=exonn($gff1);
my $exonnumber2=exonn($gff2);

foreach my $g (keys %$exonnumber1){
   print STDERR "$g\n";
   unless ($exonnumber1->{$g} == $exonnumber2->{$g}){
      print "$g\t$exonnumber1->{$g}\t$exonnumber2->{$g}\n";
   }
}



###########################
sub exonn
{
my ($gff)=@_;
my %hash;
foreach my $g (keys %$gff){
   my @line=split("\n",$gff->{$g});
   my $mrna=shift @line;
   my $exon=@line;
   #print "$g\t$exon\n";
   $hash{$g}=$exon;
}
return \%hash;
}

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
 
