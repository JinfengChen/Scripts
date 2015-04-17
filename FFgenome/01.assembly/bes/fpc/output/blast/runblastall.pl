### cut the fasta into subfiles and run blastall by qsub-sge.pl
use strict;
use warnings;
use File::Basename qw(basename dirname); 
use Getopt::Long;

my ($infile,$database,$help);
GetOptions(
         "infile:s" => \$infile,
         "database:s" => \$database,
         "help:s"    => \$help
);
die "Usage: perl runblastall.pl -i infile -d database > log" if ($help);
die "Usage: perl runblastall.pl -i infile -d database > log" unless (-f $infile);

my $filename=basename($infile);
print "$infile\t$filename\n";

my $outdir=".";
my $blastall="/share/raid1/genome/bin/blastall";
my $fastadeal="/share/raid12/chenjinfeng/tools/bin/fastaDeal.pl";
my $qsub="/share/raid12/chenjinfeng/tools/bin/qsub-sge.pl";
## cut file and push file name to array

`$fastadeal -cutf 30 $infile -outdir $outdir`;
my @subfiles=glob("./$filename.cut/*.*");

## write shell file 
my $blastshell="$filename".".sh";
open OUT, ">$blastshell" or die "can not open my shell out";
foreach (@subfiles){
    print "$_\n";
    print OUT "$blastall -p blastn -i $_ -d $database -o $_.blastm8 -e 1e-10 -m 8 > $_.log 2> $_.log2\n";
}
close OUT;

## run shell by qsub-sge.pl
`$qsub --reqsub $blastshell`;
`cat $outdir/$filename.cut/*.blastm8 > $outdir/$filename.blastm8`;
`rm -R $outdir/$filename.cut`;


