#!/usr/bin/perl
=header
The scripts is designed to run tophat to map transcriptome read to reference genome.
You can map paired and unpaired read in one run using -1,-2,-read or in two seperate run using -1,-2 and -read.
-ref:   reference sequence
-length: reference length
-gff:   gff file of gene annotation
-read:  fastq files of unpaired solexa reads seperate by ",", SRR034638.fastq,SRR034639.fastq
-1:     paired read with -2, SRR034638_1.fastq,SRR034639_1.fastq 
-2:     paired read with -1, SRR034638_2.fastq,SRR034639_2.fastq
-project: project name that used for result file 
=cut

use strict;
use warnings;
use FindBin qw ($Bin);
use File::Basename qw(basename dirname); 
use Getopt::Long;

my %opt;
GetOptions(\%opt,"ref:s","length:s","gff:s","1:s","2:s","read:s","project:s","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: perl $0 -ref ../output/OBa.all -length reflen.txt -gff FF.Gene.gff -read unpaired.fastq -1 SRR042638_1.fq -2 SRR042638_2.fq -p OBa\n";
   exit();
}

my $tophat="/share/raid1/genome/bin/tophat";
my $SAMtool="/share/raid12/chenjinfeng/tools/samtools-0.1.7_x86_64-linux/samtools";
my $bamToBed ="/share/raid12/chenjinfeng/tools/BEDTools/bin/bamToBed";
my $qsub="/share/raid12/chenjinfeng/tools/bin/qsub-sge.pl";

my $tophatshell="$opt{project}.tophat.sh";
my $tophatout="./$opt{project}"."_tophat2";
if (exists $opt{1} and exists $opt{2}){  ## paired reads
   open OUT, ">$tophatshell" or die "$!";
       print OUT "$tophat -r 50 -o $tophatout -G $opt{gff} -p 4 -g 1 --no-novel-juncs --solexa1.3-quals $opt{ref} $opt{1} $opt{2} > $opt{project}.log 2> $opt{project}.log2\n";
       #print OUT "$tophat -r 50 -o $tophatout -G $opt{gff} -p 4 -g 1 --solexa1.3-quals $opt{ref} $opt{1} $opt{2} > $opt{project}.log 2> $opt{project}.log2\n";
   close OUT;
   system ("$qsub --resource vf=0.6G $tophatshell");
}elsif(exists $opt{read}){  ## single end reads
   open OUT, ">$tophatshell" or die "$!";
       print OUT "$tophat -r 50 -o $tophatout -G $opt{gff} -p 4 -g 1 --solexa-quals $opt{ref} $opt{read} > $opt{project}.log 2> $opt{project}.log2\n";
   close OUT;
   system ("$qsub --resource vf=0.6G $tophatshell");

}

system ("$SAMtool view -bt $opt{length} $tophatout/accepted_hits.sam > $tophatout/accepted_hits.bam");
system ("$bamToBed -i $tophatout/accepted_hits.bam > $tophatout/accepted_hits.bed");
system ("rm -R $tophatshell*");


