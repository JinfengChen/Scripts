#!/usr/bin/perl
=header
The scripts is designed to run bowtie to map bisulphite sequencing read to reference genome.
You can map paired and unpaired read in one run using -1,-2,-read or in two seperate run using -1,-2 and -read.
-ref:   reference sequence
-read:  fastq files of unpaired solexa reads seperate by ",", SRR034638.fastq,SRR034639.fastq
-1:     paired read with -2, SRR034638_1.fastq,SRR034639_1.fastq 
-2:     paired read with -1, SRR034638_2.fastq,SRR034639_2.fastq
-project: project name that used for result file 
=cut

use Getopt::Long;
my %opt;
GetOptions(\%opt,"ref:s","1:s","2:s","read:s","project:s","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: perl $0 -ref /home/jfchen/epigenome/peak/input/tigr6.1/all.fa -read unpaired.fastq -1 SRR042638_1.fastq -2 SRR042638_2.fastq -project BSmethyl\n";
   exit();
}

my $bowtie="/home/jfchen/software/bowtie-0.12.5/bowtie";
my $SAMtool="/home/jfchen/software/samtools-0.1.7_x86_64-linux/samtools";
my $bamToBed ="/home/jfchen/software/BEDTools/bin/bamToBed";

## Making Bisuphite genome
#system ("mkdir bsgenome");
#my $Tgenome="./bsgenome/$opt{ref}.Tgenome.fa";
#my $Agenome="./bsgenome/$opt{ref}.Agenome.fa";
my $Tgenome="./bsgenome/Tgenome.fa";
my $Agenome="./bsgenome/Agenome.fa";
#fastaC2T($opt{ref},$Tgenome);
#fastaG2A($opt{ref},$Agenome);

## Making Bisuphite reads
system ("mkdir bsreads");
if (exists $opt{1} and exists $opt{2}){  ## paired reads
    my @read1=split(",",$opt{1});
    my @read2=split(",",$opt{2});
    my (@read1T,@read1A);
    my (@read2T,@read2A);
    for(my $i=0;$i<@read1;$i++){
        my $Tread1="./bsreads/$read1[$i].T.fastq";      
        my $Tread2="./bsreads/$read2[$i].T.fastq";
        fastqC2T($read1[$i],$Tread1);
        fastqC2T($read2[$i],$Tread2);
        push (@read1T,$Tread1);        
        push (@read2T,$Tread2);

        my $Aread1="./bsreads/$read1[$i].A.fastq";
        my $Aread2="./bsreads/$read2[$i].A.fastq"; 
        fastqG2A($read1[$i],$Aread1);
        fastqG2A($read2[$i],$Aread2);
        push (@read1A,$Aread1);
        push (@read2A,$Aread2);
    }
    my $read1T=join(",",@read1T);
    my $read1A=join(",",@read1A);
    my $read2T=join(",",@read2T);
    my $read2A=join(",",@read2A);
    ## Run bowtie
    print "Running bowtie...\n";
    system ("$bowtie -v 2 -S $Tgenome -1 $read1T -2 $read2A  $opt{project}.1T2A.T.sam > $opt{project}.1T2A.T.log");
    system ("$bowtie -v 2 -S $Agenome -1 $read1T -2 $read2A  $opt{project}.1T2A.A.sam > $opt{project}.1T2A.A.log");
    
    system ("$bowtie -v 2 -S $Tgenome -1 $read1A -2 $read2T  $opt{project}.1A2T.T.sam > $opt{project}.1A2T.T.log");
    system ("$bowtie -v 2 -S $Agenome -1 $read1A -2 $read2T  $opt{project}.1A2T.A.sam > $opt{project}.1A2T.A.log");
    print "Bowtie done!\n";
    
    ## Format changing
    print "Format changing...\n";
    system ("$SAMtool view -bS -o $opt{project}.1T2A.T.raw.bam $opt{project}.1T2A.T.sam"); 
    system ("$SAMtool view -bS -o $opt{project}.1T2A.A.raw.bam $opt{project}.1T2A.A.sam");
    system ("$SAMtool view -bS -o $opt{project}.1A2T.T.raw.bam $opt{project}.1A2T.T.sam");
    system ("$SAMtool view -bS -o $opt{project}.1A2T.T.raw.bam $opt{project}.1A2T.A.sam");
    ## remove potential PCR duplicates
    system ("$SAMtool rmdup $opt{project}.1T2A.T.raw.bam $opt{project}.1T2A.T.bam");
    system ("$SAMtool rmdup $opt{project}.1T2A.A.raw.bam $opt{project}.1T2A.A.bam");
    system ("$SAMtool rmdup $opt{project}.1A2T.T.raw.bam $opt{project}.1A2T.T.bam");
    system ("$SAMtool rmdup $opt{project}.1A2T.A.raw.bam $opt{project}.1A2T.A.bam");
    ## bam2Bed
    system ("$bamToBed -i $opt{project}.1T2A.T.bam > $opt{project}.1T2A.T.bed");
    system ("$bamToBed -i $opt{project}.1T2A.A.bam > $opt{project}.1T2A.A.bed");
    system ("$bamToBed -i $opt{project}.1A2T.T.bam > $opt{project}.1A2T.T.bed");
    system ("$bamToBed -i $opt{project}.1A2T.A.bam > $opt{project}.1A2T.A.bed");
    print "Formating done!\n";
}
if(exists $opt{read}){ ## unpaired reads
    my @read=split(",",$opt{read});
    my (@readT,@readA);
    for(my $i=0;$i<@read;$i++){  
       my $Tread="./bsreads/$read[$i].T.fastq";
       my $Aread="./bsreads/$read[$i].A.fastq";
       #fastqC2T($read[$i],$Tread);
       push (@readT,$Tread);
       #fastqG2A($read[$i],$Aread);
       push (@readA,$Aread);
    }
    my $readT=join(",",@readT);
    my $readA=join(",",@readA);    
    ## Run bowtie
    print "Running bowtie...\n";
    system ("$bowtie -v 3 -S $Tgenome $readT $opt{project}.T.T.sam > $opt{project}.T.T.log");
    system ("$bowtie -v 3 -S $Agenome $readT $opt{project}.T.A.sam > $opt{project}.T.A.log");
    system ("$bowtie -v 3 -S $Tgenome $readA $opt{project}.A.T.sam > $opt{project}.A.T.log");
    system ("$bowtie -v 3 -S $Agenome $readA $opt{project}.A.A.sam > $opt{project}.A.A.log");
    print "Bowtie done!\n";
    ## Format changing
    print "Format changing...\n";
    system ("$SAMtool view -bS -o $opt{project}.T.T.raw.bam $opt{project}.T.T.sam");
    system ("$SAMtool view -bS -o $opt{project}.T.A.raw.bam $opt{project}.T.A.sam");
    system ("$SAMtool view -bS -o $opt{project}.A.T.raw.bam $opt{project}.A.T.sam");
    system ("$SAMtool view -bS -o $opt{project}.A.A.raw.bam $opt{project}.A.A.sam");
    ## remove potential PCR duplicates
    system ("$SAMtool rmdup -s $opt{project}.T.T.raw.bam $opt{project}.T.T.bam");
    system ("$SAMtool rmdup -s $opt{project}.T.A.raw.bam $opt{project}.T.A.bam");
    system ("$SAMtool rmdup -s $opt{project}.A.T.raw.bam $opt{project}.A.T.bam");
    system ("$SAMtool rmdup -s $opt{project}.A.A.raw.bam $opt{project}.A.A.bam");
    ## bam2Bed
    system ("$bamToBed -i $opt{project}.T.T.bam > $opt{project}.T.T.bed");
    system ("$bamToBed -i $opt{project}.T.A.bam > $opt{project}.T.A.bed");
    system ("$bamToBed -i $opt{project}.A.T.bam > $opt{project}.A.T.bed");
    system ("$bamToBed -i $opt{project}.A.A.bam > $opt{project}.A.A.bed");
    print "Formating done!\n";
}





sub fastqC2T
{
my ($infile,$outfile)=@_;
open IN, "$infile" or die "$!";
open OUT, ">$outfile" or die "$!";
while (<IN>){
    if ($_=~/^@/){
       print OUT "$_";
       my $seq=<IN>;
       $seq=~tr/Cc/Tt/;
       print OUT "$seq";
       my $head1=<IN>;
       my $qual =<IN>;
       print OUT "$head1$qual";
    }
}
close OUT;
close IN;
}

sub fastqG2A
    {
    my ($infile,$outfile)=@_;
    open IN, "$infile" or die "$!";
    open OUT, ">$outfile" or die "$!";
    while (<IN>){
      if ($_=~/^@/){
        print OUT "$_";
        my $seq=<IN>;
        $seq=~tr/Gg/Aa/;
        print OUT "$seq";
        my $head1=<IN>;
        my $qual =<IN>;
        print OUT "$head1$qual";
      }
    }
    close OUT;
    close IN;
}



sub fastaC2T
{
my ($infile,$outfile)=@_;
$/=">";
open IN, "$infile" or die "$!";
open OUT, ">$outfile" or die "$!";
while (<IN>){
   next if (length $_ <= 2);
   my @unit=split("\n",$_);
   my $head=shift @unit;
   my $seq =join("\n",@unit);
   $seq=~tr/Cc/Tt/;
   $seq=~s/>//g;
   print "$head\n";
   print OUT ">$head\n$seq\n";
}
close OUT;
close IN;
}
sub fastaG2A
{
     my ($infile,$outfile)=@_;
     $/=">";
     open IN, "$infile" or die "$!";
     open OUT, ">$outfile" or die "$!";
     while (<IN>){
         next if (length $_ <= 2);
         my @unit=split("\n",$_);
         my $head=shift @unit;
         my $seq =join("\n",@unit);
         $seq=~tr/Gg/Aa/;
         $seq=~s/>//g;
         print OUT ">$head\n$seq\n";
     }
       close OUT;
       close IN;
}
