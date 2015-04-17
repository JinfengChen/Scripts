#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"gff1:s","gff2:s","orth:s","syntid:s","help");


my $help=<<USAGE;
Some ortholog genes are tend to be a duplication of original ortholog and present on other chromosome.
This scripts is to find these ortholog pairs.
perl $0 --gff1 -gff2 --orth --syntid
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $refgff1=parseGFF($opt{gff1});
my $refgff2=parseGFF($opt{gff2});
my $refid  =id($opt{syntid}); ## family id
my $reforth=ortholog($opt{orth});

foreach my $gene1 (keys %$reforth){
    my $flag=0;
    ##geenid chr startposition
    my $ref="$gene1\t$refgff2->{$gene1}->[0]\t$refgff2->{$gene1}->[1]";
    my @qry;
    foreach my $refgene (@{$reforth->{$gene1}}){
        my $gene2=$refgene->[0];
        my $identity=$refgene->[1];
        push (@qry,"$gene2\t$refgff1->{$gene2}->[0]\t$refgff1->{$gene2}->[1]\t$identity");
        if ($refgff2->{$gene1}->[0] ne $refgff1->{$gene2}->[0]){
           $flag=1;
        }
    }
    if ($flag == 1){
       my $line=join("\t",@qry);
       print "$ref\t$refid->{$gene1}\t$line\n";
    }
}


#############################

sub id
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $unit[0]=$1 if ($unit[0]=~/(.*)\_\w+$/);
    $hash{$unit[0]}=$unit[1];
}
close IN;
return \%hash;
}

sub ortholog
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    ## gene1 cdsidentity
    push (@{$hash{$unit[1]}},[$unit[0],$unit[3]]);
}
close IN;
return \%hash;
}



sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;   ##Scaffold
my $id;    ##ID for element
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*?);/){
            $id=$1;
        }
        $hash{$id}=[$seq,$unit[3]];
    }

}
close IN;
return \%hash;
}

