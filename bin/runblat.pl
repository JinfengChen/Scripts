#!/usr/bin/perl
use strict;
use warnings;
use FindBin qw ($Bin);
use File::Basename qw(basename dirname); 
use Getopt::Long;

my ($infile,$evalue,$database,$blast,$help);
GetOptions(
         "infile:s" => \$infile,
         "database:s" => \$database,
         "help:s"    => \$help
);

unless (-e $infile and -e $database){
   print "input file speficied!\n";
   exit();
}

my $filename=basename($infile);
print "$infile\t$filename\n";

my $fastadeal="/rhome/cjinfeng/software/bin/fastaDeal.pl";
my $qsub="/rhome/cjinfeng/software/bin/qsub-pbs.pl";
my $bestalign="/rhome/cjinfeng/software/bin/bestAlign.pl";

my $outdir="./";
my $fdshell="$filename".".cut".".sh";
open OUT, ">$fdshell" or die "$!";
     print OUT "$fastadeal -cutf 60 $infile -outdir $outdir\n";
close OUT;
`perl $qsub $fdshell`;


my @subfiles=glob("$outdir/$filename.cut/*.*");

## write shell file 
my $blastshell="$filename".".sh";
open OUT, ">$blastshell" or die "can not open my shell out";
foreach (@subfiles){
    print "$_\n";
    print OUT "blat -noHead $database $_ $_.psl > $_.log 2> $_.log2\n";
    print OUT "perl $bestalign $_.psl > $_.best.psl\n";
}
close OUT;

## run shell by qsub-sge.pl
`perl $qsub --lines 2 $blastshell`;
#`cat $outdir/$filename.cut/*.psl > $outdir/$filename.psl`;
`cat $outdir/$filename.cut/*.best.psl > $outdir/$filename.best.psl`;
#rm -R $outdir/$filename.cut`;


