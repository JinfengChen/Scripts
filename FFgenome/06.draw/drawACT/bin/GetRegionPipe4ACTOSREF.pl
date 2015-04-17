#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"block:s","help");


my $help=<<USAGE;
Before run this pipe, read blockinf and change some region that have reverse direction,OS_chr01_24189665_21235376
Also, delete some regions that preduced by small align
--block: blockinf
Run: perl $0 --block testblock
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $chr;
my $outdir="../output/4ACT";
`mkdir $outdir` unless (-f $outdir);
open IN, "$opt{block}" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    if ($_=~/^(chr\d+)/){
       $chr=$1;
       $outdir1="$outdir/$chr";
       `mkdir $outdir1` unless (-f $outdir1);
       next;
    }
    my @unit=split("\t",$_);
    #print "$unit[0]\t$unit[1]\n";
    unless ($unit[0]=~/\w+\_/ and $unit[1]=~/\w+\_/){
       next;
    }
    `perl GetRegionData.pl --qryhead $unit[0] --qrygff3 ../input/9311.glean.final.gff.chr --qryfasta ../input/9311.genome.final.fa.chr --qryrepeat ../input/9311.genome.final.fa.RepeatMasker.out.gff.chr --refhead $unit[1] --refgff3 ../input/tigr.all.final.gff.chr --reffasta ../input/all.con.chr -refrepeat ../input/all.con.RepeatMasker.out.gff.chr --outdir $outdir1`;
}
close IN;
