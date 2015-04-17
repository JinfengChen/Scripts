#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"transcript:s","gene:s","genome:s","project:s","help");


my $help=<<USAGE;
perl $0 --transcript --gene --genome --project
Summary genome/gene coverage based on map rate of transcript.
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

`blat $opt{genome} $opt{transcript} RNAseq2genome.500.blat.$opt{project}.psl -noHead -minIdentity=95 > log 2> log2`;
`blat $opt{gene} $opt{transcript} RNAseq2gene.500.blat.$opt{project}.psl -noHead -minIdentity=95 > log 2> log2`;

`perl /home/jfchen/FFproject/tools/bin/bestAlign.pl RNAseq2genome.500.blat.$opt{project}.psl --cutoff 0.3 > RNAseq2genome.500.blat.$opt{project}.0.3.psl`;
`perl /home/jfchen/FFproject/tools/bin/bestAlign.pl RNAseq2gene.500.blat.$opt{project}.psl --cutoff 0.3 > RNAseq2gene.500.blat.$opt{project}.0.3.psl`;

my $all=`grep ">" -c $opt{transcript}`;
my $ingene1=`cut -f10 RNAseq2gene.500.blat.$opt{project}.psl | uniq | sort | uniq | wc -l`;
my $ingenome1=`cut -f10 RNAseq2genome.500.blat.$opt{project}.psl | uniq | sort | uniq | wc -l`;
my $ingene2=`cut -f10 RNAseq2gene.500.blat.$opt{project}.0.3.psl | uniq | sort | uniq | wc -l`;
my $ingenome2=`cut -f10 RNAseq2genome.500.blat.$opt{project}.0.3.psl | uniq | sort | uniq | wc -l`;
chomp ($ingene1,$ingenome1,$ingene2,$ingenome2,$all);
#print "$all\t$ingene\t$ingenome\n";
my $generate1=$ingene1/$all;
my $generate2=$ingene2/$all;
my $genomerate1=$ingenome1/$all;
my $genomerate2=$ingenome2/$all;
print "All\t$all\n";
print "Genome Coverage\t$ingenome1($genomerate1)\t$ingenome2($genomerate2)\n";
print "Gene Coverage\t$ingene1($generate1)\t$ingene2($generate2)\n";



