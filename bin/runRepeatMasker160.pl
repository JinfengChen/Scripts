use strict;
use warnings;
use FindBin qw($Bin);
use File::Basename qw(basename dirname);


die "Usage: perl runrepeatmasker.pl infile lib > log" if (@ARGV < 1);
my $infile=$ARGV[0];
my $lib ||= $ARGV[1];
my $filename=basename($infile);
print "$infile\t$filename\n";

my $outdir="./";
my $common="/home/jfchen/159/FFproject/tools/bin";
my $repeat2gff="$common/repeat_to_gff.pl";
my $fastadeal="$common/fastaDeal.pl";
my $repeatmasker="RepeatMasker";
my $statTE="$common/stat_TE.pl";
my $qsub="$common/multi-process.pl";
## cut file and push file name to array
my $fdshell="$filename".".cut".".sh";
open OUT, ">$fdshell" or die "$!";
     print OUT "$fastadeal -cutf 20 $infile -outdir $outdir\n";
close OUT;
`$qsub $fdshell`;

my @subfiles=glob("$outdir/$filename.cut/*.*");

## write shell file 
my $repeatshell="$filename".".sh";
open OUT, ">$repeatshell" or die "can not open my shell out";
foreach (@subfiles){
    print "$_\n";
    unless (-f $lib){
       print OUT "$repeatmasker -species rice -xsmall -nolow -no_is -qq -norna -engine wublast $_ > $_.log 2> $_.log2\n";
    }else{
       print OUT "$repeatmasker -lib $lib -nolow -xsmall -no_is -norna -engine wublast $_ > $_.log 2> $_.log2\n";
    }
}
close OUT;

## run shell by qsub-sge.pl
`$qsub --cpu 1 $repeatshell`;
`cat $outdir/$filename.cut/*.out > $outdir/$filename.RepeatMasker.out`;
`cat $outdir/$filename.cut/*.masked > $outdir/$filename.RepeatMasker.masked`;
`cat $outdir/$filename.cut/*.tbl > $outdir/$filename.RepeatMasker.tbl`;
`cat $outdir/$filename.cut/*.cat > $outdir/$filename.RepeatMasker.cat`;
`perl $repeat2gff $outdir/$filename.RepeatMasker.out`;
`perl $statTE --repeat $outdir/$filename.RepeatMasker.out --rank all > $outdir/$filename.RepeatMasker.out.stat.all`;
`perl $statTE --repeat $outdir/$filename.RepeatMasker.out --rank type > $outdir/$filename.RepeatMasker.out.stat.type`;
`perl $statTE --repeat $outdir/$filename.RepeatMasker.out --rank subtype > $outdir/$filename.RepeatMasker.out.stat.subtype`;
`perl $statTE --repeat $outdir/$filename.RepeatMasker.out --rank family > $outdir/$filename.RepeatMasker.out.stat.family`;
`rm -R $outdir/$filename.cut`;

