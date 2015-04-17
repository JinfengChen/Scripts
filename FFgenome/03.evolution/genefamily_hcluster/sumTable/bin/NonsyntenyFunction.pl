#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"gff:s","id:s","iprscan:s","help");


my $help=<<USAGE;
perl $0 --gff --id
count gene number for each annotation function.
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $refgff=parseGFF($opt{gff});
my ($refipr,$ipr)=readiprscan($opt{iprscan});
rankanno($opt{id},$refgff,$refipr,$ipr);
#######

sub rankanno
{
my ($file,$gff,$refipr,$ipr)=@_;
my %hash;
my %pfam;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $unit[0]=$1 if ($unit[0]=~/(.*)_\w+/);
    $hash{$gff->{$unit[0]}}++;
    $pfam{$refipr->{$unit[0]}}++;
}
close IN;
#foreach my $anno (sort {$hash{$b} <=> $hash{$a}} keys %hash){
    #print "$anno\t$hash{$anno}\n";
#}
foreach my $pfam (sort {$pfam{$b} <=> $pfam{$a}} keys %pfam){
    print "$pfam\t$ipr->{$pfam}\t$pfam{$pfam}\n";
}
}

#################
sub readiprscan{
my ($file)=@_;
my %hash;
my %id;
my %ipr;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    next unless ($unit[3] eq "HMMPfam");
    next if (exists $id{$unit[0]});
    next if ($unit[11] eq "NULL");
    $ipr{$unit[11]}=$unit[12];
    $hash{$unit[0]}= $unit[11];
}
close IN;
return (\%hash,\%ipr);
}


sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $id;    ##ID for element
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*?);Name=(.*?);/){
            $id=$1;
            $hash{$id}=$2;
        }
    }

}
close IN;

return \%hash;
}
 
