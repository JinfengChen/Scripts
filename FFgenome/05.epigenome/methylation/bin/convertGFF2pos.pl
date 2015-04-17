#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"gff:s","project:s","help");


my $help=<<USAGE;

Convert gff format file to bed format.
Example:
Scaffold000021  4057231 4057449 gene    OBR_GLEAN_10015563      +

Run: perl $0 -gff OBa.all.gff -project FF_gene_bed
-gff: GFF file 
-project: prefix for output file

USAGE

if (keys %opt < 1){
   print "$help\n";
   exit();
}

`mkdir $opt{project}` unless (-f $opt{project});
 
my $refgff=parseGFF($opt{gff});
foreach(keys %$refgff){
    my @array=split("\n",$refgff->{$_});
    my $gene=$_;
    my @pos;
    my $chr;
    my $strand;
    foreach(@array){
       my @unit=split("\t",$_);
       if ($unit[2] eq "CDS"){
          push (@pos,$unit[3]);
          push (@pos,$unit[4]);
       }elsif($unit[2] eq "mRNA"){
          $chr=$unit[0];
          $strand=$unit[6];
       }
    }
    my @temp=sort {$a <=> $b} @pos;
    my $start=shift @temp;
    my $end=pop @temp;
    open OUT, ">>$opt{project}/$chr.gene.bed" or die "$!";
    print OUT "$chr\t$start\t$end\tgene\t$gene\t$strand\n";
    close OUT;
}


##########################################
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


