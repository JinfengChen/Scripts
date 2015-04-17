#!/usr/bin/perl
### cut the fasta into subfiles and run blastall by qsub-sge.pl
### be sure: cp the script into your project bin directory or where you run it
use strict;
use warnings;
use FindBin qw ($Bin);
use File::Basename qw(basename dirname); 
use Getopt::Long;

my ($blastn, $infile,$database,$help);
GetOptions(
         "infile:s" => \$infile,
         "blastn:s"   => \$blastn,
         "database:s" => \$database,
         "help:s"    => \$help
);
die "Usage: perl runblastall.pl -b blastn -i infile -d database > log" if ($help);
die "Usage: perl runblastall.pl -b blastn -i infile -d database > log" unless (-f $infile);

my $filename=basename($infile);
print "$infile\t$filename\n";

my $outdir= `pwd`;
chomp $outdir;
my $blastall="/usr/bin/blastall";
my $fastadeal="/rhome/cjinfeng/software/bin/fastaDeal.pl";
my $qsub="/rhome/cjinfeng/software/bin/qsub-pbs.pl";
## cut file and push file name to array
#`$fastadeal -cutf 60 $infile -outdir $outdir`;


my $fdshell="$filename".".cut".".sh";
open OUT, ">$fdshell" or die "$!";
     print OUT "perl $fastadeal -cutf 40 -outdir $outdir $infile\n";
close OUT;
`perl $qsub $fdshell`;


my @subfiles=glob("$outdir/$filename.cut/*.*");

## write shell file 
my $blastshell="$filename".".sh";
open OUT, ">$blastshell" or die "can not open my shell out";
foreach (@subfiles){
    #print "$_\n";
    print OUT "$blastall -p $blastn -i $_ -d $database -o $_.blastm8 -e 1e-5 -m 8 > $_.log 2> $_.log2\n";
}
close OUT;

## run shell by qsub-sge.pl
`perl $qsub --resource nodes=1:ppn=1,mem=3G,walltime=100:00:00 $blastshell`;
`cat $outdir/$filename.cut/*.blastm8 > $outdir/$filename.blastm8`;
`rm -R $outdir/$filename.cut`;


