#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"iprscan:s","hcluster:s","tandem:s","ortholog:s","gff:s","project:s","help");


my $help=<<USAGE;
Annotate Pfam family and gene family id for each tandem cluster.
Run: perl $0 --iprscan ../input/representative_orf.fa.nr.pep.iprscan --tandem OS.tandem_repeat_gene.txt --hcluster ../input/all_vs_all.blast.m8.solar.forHC.hcluster --ortholog ../input/IRGSP.ortholog.table --project IRGSP > log 2> log2 &
--iprscan: iprscan result file
--tandem: 
--hcluster
--ortholog:
--gff:
--project: suffix of gene in ortholog_paralog_analysis, OBR_GLEAN_10025455_OBRACH (OBRACH)
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $refiprscan=iprscan($opt{iprscan});
my $reffamily=family($opt{hcluster});
my $reforth=ortholog($opt{ortholog});
my $refstart=genestart($opt{gff});

& tandem($opt{tandem},$refiprscan,$reffamily,$reforth,$refstart);

#############################
sub iprscan
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ =~/^$/);
    my @unit=split("\t",$_);
    next if (exists $hash{$unit[0]});
    next unless ($unit[3]=~/HMMPfam/);
    $hash{$unit[0]}=[$unit[4],$unit[5]];
}
close IN;
return \%hash;
}

sub family
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ =~/^$/);
    my @unit=split("\t",$_);
    my @gene=split(",",$unit[6]);
    foreach my $gene (@gene){
       if ($gene=~/(.*)\_$opt{project}/){
          $hash{$1}=$unit[0];
       }
    }
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
    next if ($_ =~/^$/);
    my @unit=split("\t",$_);
    if ($unit[0]=~/(.*)\_$opt{project}/){
        $hash{$1}=1;
    }
}
close IN;
return \%hash;
}

sub genestart
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
        if ($unit[8]=~/ID=(.*\.\d+);/ or $unit[8] =~/ID=(.*\.\d+)/){
            $id=$1;
        }
        $hash{$id}=$unit[3]; 
    }
}
close IN;
return \%hash;
}


sub tandem
{
my ($file,$refiprscan,$reffamily,$reforth,$refstart)=@_;
my $counter;
open OUT, ">$opt{project}.tandem.anno";
print OUT "TandemID\tGeneNumber\tSyntenyGeneNumber\tPosition\tFamilyID\tPfamAnnotation\tGeneIDs\tFamilyIDs\tPfamIDs\tPfamAnnotations\n";
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    $counter++;
    my $line=$_;
    my @unit=split(" ",$_);
    my $tandem=join(",",@unit);
    my @family;
    my @pfam;
    my @pfamid;
    my @start;
    my $genen=0;
    my $orthn=0;
    foreach my $gene (@unit){
        $genen++;
        $orthn++ if (exists $reforth->{$gene});
        my $fam= $reffamily->{$gene} ? $reffamily->{$gene} : "NA"; 
        my $pfamid= $refiprscan->{$gene}->[0] ? $refiprscan->{$gene}->[0] : "NA";
        my $pfam= $refiprscan->{$gene}->[1] ? $refiprscan->{$gene}->[1] : "NA";
        my $start= $refstart->{$gene} ? $refstart->{$gene} : "NA";
        push (@family,$fam);
        push (@pfam,$pfam);
        push (@pfamid,$pfamid);
        push (@start,$start);
    }
    my $fam1=join(",",@family);
    my $pfam1=join(",",@pfam);
    my $pfamid1=join(",",@pfamid);
    print OUT "$counter\t$genen\t$orthn\t$start[0]\t$family[0]\t$pfam[0]\t$tandem\t$fam1\t$pfamid1\t$pfam1\n";
}
close IN;
close OUT;
}
