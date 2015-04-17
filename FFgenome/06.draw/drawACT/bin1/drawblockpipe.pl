#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"block:s","help");


my $help=<<USAGE;
--block: blockinf
Run: perl $0 --block testblock
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $chr;
my $outdir="../output";
my $pdf="../output/pdf";
`mkdir $pdf` unless (-f $pdf);
open IN, "$opt{block}" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    if ($_=~/^(chr\d+)/){
       $chr=$1;
       $outdir="../output/$chr";
       $pdf="../output/pdf/$chr";
       `mkdir $outdir` unless (-f $outdir);
       `mkdir $pdf` unless (-f $pdf);
       next;
    }
    my @unit=split("\t",$_);
    print "$unit[0]\t$unit[1]\n";
    unless ($unit[0]=~/\w+\_/ and $unit[1]=~/\w+\_/){
       next;
    }
    `perl GetRegionData.pl --qryhead $unit[0] --qrygff3 ../data/Ob.gff3.chr --qryfasta ../data/Ob.fasta.chr --qryrepeat ../data/OBa.all.fa.RepeatMasker.out.gff.chr --refhead $unit[1] --refgff3 ../data/Os.gff3.chr --reffasta ../data/Os.fasta.chr --refrepeat ../data/IRGSP.build5.RepeatMasker.out.gff.chr --outdir $outdir --qRNAseq ../data/Ob.Shoot.bed.chr --qmC ../data/FF.methylation.bed.chr --qH3K4 ../data/OB.H3K4.quarter.bed.chr --rRNAseq ../data/IRGSPShoot.bed.chr --rmC ../data/IRGSP.methylation.bed.chr --rH3K4 ../data/OS.H3K4.bed.chr`;
    `perl drawACTgraph.pl -query $outdir/$unit[1]/$unit[1].embl -target $outdir/$unit[1]/$unit[0].embl -act $outdir/$unit[1]/$unit[1]VS$unit[0]4ACT -desc $unit[1]vs$unit[0] -qRNAseq $outdir/$unit[1]/$unit[1].RNAseq.bed -qmC $outdir/$unit[1]/$unit[1].mC.bed -qH3K4 $outdir/$unit[1]/$unit[1].H3K4.bed -tRNAseq $outdir/$unit[1]/$unit[0].RNAseq.bed -tmC $outdir/$unit[1]/$unit[0].mC.bed -tH3K4 $outdir/$unit[1]/$unit[0].H3K4.bed > draw.log 2> draw.log2`;
    `mv $unit[1]* $pdf`;
}
close IN;
