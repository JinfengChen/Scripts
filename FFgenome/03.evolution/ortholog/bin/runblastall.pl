#!/usr/bin/perl
### cut the fasta into subfiles and run blastall by qsub-sge.pl
### be sure: cp the script into your project bin directory
use strict;
use warnings;
use FindBin qw ($Bin);
use File::Basename qw(basename dirname); 
use Getopt::Long;

my ($infile,$format,$database,$blast,$help);
GetOptions(
         "program:s"  => \$blast,
	 "format:s"  => \$format, ### T if protein, F if nucleotide
         "infile:s" => \$infile,
         "database:s" => \$database,
         "help:s"    => \$help
);

$blast ||="blastn";

die "Usage: perl runblastall.pl -p blastp -f T -i infile -d database > log" if ($help);
die "Usage: perl runblastall.pl -p blastp -f T -i infile -d database > log" unless (-f $infile);

my $filename=basename($infile);
my $dbname=basename($database);
print "$infile\t$filename\n";

my $outdir="../output";
my $blastall="blastall";
my $formatdb="formatdb";
my $fastadeal="/home/jfchen/FFproject/tools/bin/fastaDeal.pl";
my $qsub="/home/jfchen/FFproject/tools/bin/multi-process.pl";
my $blastparser="/home/jfchen/FFproject/tools/bin/blast_parser.pl";
## cut file and push file name to array
#`$fastadeal -cutf 60 $infile -outdir $outdir`;


my $fdshell="$filename".".cut".".sh";
open OUT, ">$fdshell" or die "$!";
     print OUT "$fastadeal -cutf 6 $infile -outdir $outdir\n";
close OUT;
#`$qsub --resource vf=0.2G $fdshell`;
`$qsub --cpu 6 $fdshell`;

my @subfiles=glob("$outdir/$filename.cut/*.*");
#system("$formatdb -p $format -i $database");

## write shell file 
my $blastshell="$filename".".sh";
open OUT, ">$blastshell" or die "can not open my shell out";
foreach (@subfiles){
    print "$_\n";
    print OUT "$blastall -p $blast -F F -i $_ -d $database -o $_.blast -e 1e-10 > $_.log 2> $_.log2\n";
}
close OUT;

my $blastout=$filename."_".$dbname.".blast";
my $parseout=$filename."_".$dbname.".blasttable";
my $inparanoid=$filename."-".$dbname;
## run shell by qsub-sge.pl
`$qsub --cpu 6 $blastshell`;
`cat $outdir/$filename.cut/*.blast > $outdir/$blastout`;
`$blastparser $outdir/$blastout > $outdir/$parseout`;
#`$table2inparanoid $outdir/$parseout > $inparanoid`;
`rm -R $outdir/$filename.cut`;
`rm $outdir/$blastout`;


