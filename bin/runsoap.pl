#!/usr/bin/perl
### cut the fasta into subfiles and run soap by qsub-sge.pl
use strict;
use warnings;
use FindBin qw ($Bin);
use File::Basename qw(basename dirname); 
use Getopt::Long;

my ($infile,$database,$help);
GetOptions(
         "infile:s" => \$infile,
         "database:s" => \$database,
         "help:s"    => \$help
);
die "Usage: perl runsoap.pl -i infile -d database.index > log &" if ($help);
die "Usage: perl runsoap.pl -i infile -d database.index > log &" unless (-f $infile);

my $filename=basename($infile);
#print "$infile\t$filename\n";

my $outdir="$Bin/../output";
#print "$outdir\n";

my $soap="/share/raid12/chenjinfeng/tools/soapall/soap2.17release/soap";
my $fastadeal="/share/raid12/chenjinfeng/tools/bin/fastaDeal.pl";
my $qsub="/share/raid12/chenjinfeng/tools/bin/qsub-sge.pl";
## cut file and push file name to array
my $fdshell="$filename".".cut".".sh";
open OUT, ">$fdshell" or die "$!";
     print OUT "$fastadeal -cutf 60 $infile -outdir $outdir\n";
close OUT;
`$qsub $fdshell`;

my @subfiles=glob("$outdir/$filename.cut/*.*");

## write shell file 
my $soapshell="$filename".".sh";
open OUT, ">$soapshell" or die "can not open my shell out";
foreach (@subfiles){
    print "$_\n";
    print OUT "$soap -a $_ -D $database -o $_.soapout > $_.log 2> $_.log2\n";
}
close OUT;

## run shell by qsub-sge.pl
`$qsub $soapshell`;
`cat $outdir/$filename.cut/*.soapout > $outdir/$filename.soapout`;
`rm -R $outdir/$filename.cut`;

