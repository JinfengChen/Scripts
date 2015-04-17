#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"block:s","help");


my $help=<<USAGE;
Summary collinearity from block file
perl $0 --block
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

$chr= $1 if ($opt{block}=~/(.*)\.blocks/);

my $refbed=readbed("$chr.bed");
& parseblock($opt{block},$refbed);
#########

sub readbed
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[3]}=$unit[1];
}
close IN;
return \%hash;
}



sub parseblock
{
my ($block,$bed)=@_;
my ($refer,$refn,$genen);
my (%ref,%qry,$pair);
my (%refgene,%qrygene);
open OUT, ">$block.table" or die "$!";
open OUT1, ">$block.geneinblock" or die "$!";
open OUT2, ">$block.summary";
open IN, "$block" or die "$!";
while(<IN>){
     chomp $_;
     next if ($_=~/^$/);
     my $line=$_;
     if ($line=~/## View (\d+)\: pivot (\w+)$/){
        $refn=$1;
        $refer=$2; 
     }elsif($line=~/^#/){
        next;
     }elsif($line=~/^\s+(\d+)\-\s*(\d+)\:\s+?(.*)$/){
        my $genen=$2;
        my $genes=$3;
        my @unit=split("\t",$genes);
        if ($unit[0]=~/\w+/ and $unit[1]=~/\w+/){
            $pair++;
            #print OUT "$refer\t$refn\t$genen\t$unit[0]\t$unit[1]\n";
            my @gene1=split(";",$unit[0]);
            my @gene2=split(";",$unit[1]);
            my $size1=@gene1;
            my $size2=@gene2;
            ##chr gene1 gene1pos tandemnumber gene2 gene2pos tandemnumber
            if ($size1 > 1 or $size2 > 1){
               print OUT "$refer\t$gene1[0]\t$bed->{$gene1[0]}\t$size1\t$gene2[0]\t$bed->{$gene2[0]}\t$size2\n";
            }
            foreach my $g1 (@gene1){
               $ref{$g1}=1;
            }
            foreach my $g2 (@gene2){
               $qry{$g2}=2;
            }
        }
        if($unit[0]=~/\w+/){
            my @gene1=split(";",$unit[0]);
            foreach my $g1 (@gene1){
               $refgene{$g1}=1;
               print OUT1 "$g1\n";
            }
        }
        if($unit[1]=~/\w+/){
            my @gene2=split(";",$unit[1]);
            foreach my $g2 (@gene2){
               $qrygene{$g2}=1;
            }
        }
     } 
}
close IN;
close OUT;
close OUT1;
my $refgenen=keys %ref;
my $qrygenen=keys %qry;
my $refgenet=keys %refgene;
my $qrygenet=keys %qrygene;
#print OUT2 "Total Pair: $pair\n";
#print OUT2 "Reference Gene Number: $refgenen\n";
#print OUT2 "Target Gene Number: $qrygenen\n";
#print OUT2 "Total Ref Gene Number: $refgenet\n";
#print OUT2 "Total Qry Gene Number: $qrygenet\n";
print OUT2 "Chr\t#Pair\t#RefGeneNumber\t#TarGeneNumber\t#TotalRefGeneNumber\t#TotalTarGeneNumber\n";
print OUT2 "$chr\t$pair\t$refgenen\t$qrygenen\t$refgenet\t$qrygenet\n";
close OUT2;
}## end of sub


