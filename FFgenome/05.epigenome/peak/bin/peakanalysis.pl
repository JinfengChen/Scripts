#!/usr/bin/perl
=header
The scripts is designed to run bowtie to align short read to reference genome, then call peak using MACS.
Some general analysis can be down with options.
-ref:   reference sequence
-read:  fastaq files of solexa reads seperate by ",", SRR034638.fastq,SRR034639.fastq
-tsize: tag size, solexa read length
-gsize: genome size of reference genome, tigr6 372317567, FF 260822143, IRGSP 382788128
-step: 1 Run bowtie to map read to genome, 2 Format changing, 3 Run MACS, 4 Change wig to bar
-project: project name that used for result file 
-flat:  refFlat.txt file of gene annotation infomations, optional
=cut

use Getopt::Long;

my %opt;
GetOptions(\%opt,"ref:s","read:s","step:s","project:s","tsize:s","gsize:s","flat:s","help");

if ($opt{help} or keys %opt < 1){
   print "Usage:In fastq dir, run perl $0 -ref /home/jfchen/epigenome/peak/input/tigr6.1/all.fa -read SRR034648.fastq -tsize 35 -gsize 372317567 -project H3K27 -step 1234 -flat refFlat.txt\n";
   exit();
}

my $bowtie="/home/jfchen/software/bowtie-0.12.5/bowtie";
my $macs  ="/home/jfchen/software/python/bin/macs";
my $SAMtool="/home/jfchen/software/samtools-0.1.7_x86_64-linux/samtools";
my $bamToBed ="/home/jfchen/software/BEDTools/bin/bamToBed";
my $rmdup="/home/jfchen/software/picard-tools-1.31/MarkDuplicates.jar";

## Run bowtie
if ($opt{step} =~/1/){
print "Running bowtie...\n";
system ("$bowtie -S $opt{ref} $opt{read} $opt{project}.sam");
print "Bowtie done!\n";
}

## Format changing
if ($opt{step} =~/2/){
print "Format changing...\n";
system ("$SAMtool view -bS -o $opt{project}.raw.bam $opt{project}.sam");
## remove potential PCR duplicates
system ("$SAMtool sort $opt{project}.raw.bam $opt{project}.sort");
system ("java -jar $rmdup ASSUME_SORTED=TRUE REMOVE_DUPLICATES=TRUE INPUT=$opt{project}.sort.bam OUTPUT=$opt{project}.bam METRICS_FILE=$opt{project}.dupli");
system ("$bamToBed -i $opt{project}.bam > $opt{project}.bed");
print "Formating done!\n";
}

## Run MACS
if ($opt{step} =~/3/){
print "Running MACS...\n";
system ("$macs -t $opt{project}.bed --mfold=10 --tsize $opt{tsize} --gsize $opt{gsize} --name $opt{project} --wig --space 50 > $opt{project}.log 2> $opt{project}.log2");
print "MACS done!\n";
}

if ($opt{step} =~/4/){
## Change wig file to bar.txt, which could be view by cisgenome
my $wigdir="./$opt{project}_MACS_wiggle/treat";
wig2bar($wigdir,$opt{project});
}


##### subfunctions##############################################################

sub wig2bar
{
my ($dir,$name)=@_;
print "$dir\n";
my $bar="$name.bar.txt";
open OUT, ">$bar" or die "$!";
while(glob("$dir/*.*")){
   print "$_\n";
   my $chr;
   my $file=$_;
   if ($file=~/(.*_(.*))\.wig\.gz/){
      system ("gunzip -d $file");
      $chr=$2;
      my $newfile="$1.wig";
      open IN, "$newfile" or die "$!";
         <IN>;
         <IN>;
         while (<IN>){
             next if ($_=~/^$/);
             print OUT "$chr\t$_";
         }
      close IN;
   }elsif($file=~/(.*_(.*))\.wig/){
      $chr=$2;
      print "$file\n";
      open IN, "$file" or die "$!";
         <IN>;
         <IN>;
         while (<IN>){
            next if ($_=~/^$/);
            print OUT "$chr\t$_";
         }
     close IN;
   }
}
close OUT;

}
