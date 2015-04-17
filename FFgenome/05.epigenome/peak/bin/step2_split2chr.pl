#!/usr/bin/perl
=header
The script is designed to split BED format of bowtie map result, as well as fastq reads, into subfiles according to the chromosome they are mapped. Then, run step2_sumCm.pl to summary methylation.
BED files from seperate runs of step1 should be merge into T genome BED (BSmethyl.T.bed) and A genome BED (BSmethyl.A.bed).
Value of -1,-2,-read should be same as in step1.
=cut
use Getopt::Long;
our %opt;
GetOptions(\%opt,"ref:s","1:s","2:s","read:s","project:s","help");

if ($opt{help} or keys %opt < 1){
   print "Usage: perl $0 -ref /home/jfchen/epigenome/peak/input/tigr6.1/all.fa -1 SRR042638_1.fastq -2 SRR042638_2.fastq -read SRR042639_45.fastq,SRR042639_31.fastq -project BSmethyl\n";
   exit();
}

#our (%readT);  ### globle variable, store read->chr relation for T genome
#our (%readA);  ### globle variable, store read->chr relation for A genome
our (%chr);   ### globle variable, store chr name
system ("mkdir temp_chr") unless (-e temp_chr);
### split BED
my $readT=splitBED($opt{project},".T.bed");
my $readA=splitBED($opt{project},".A.bed");

### split fastq reads
if (exists $opt{1} and exists $opt{2}){  ## paired reads
    my @read1=split(",",$opt{1});
    my @read2=split(",",$opt{2});
    splitread(\@read1,1,$readT,$readA);
    splitread(\@read2,2,$readT,$readA);
    
}elsif(exists $opt{read}){  ## unpaired reads
    my @read=split(",",$opt{read}); 
    splitread(\@read,0,$readT,$readA);
}

### run step2_sumCm.pl for each chr and merge the results
if (exists $opt{1} and exists $opt{2} and exists $opt{read}){  ## paired and unpaired reads
   foreach (keys %chr){
       my $project="./temp_chr/$_";
       my $read1="./temp_chr/"."$_"."_1".".fastq";
       my $read2="./temp_chr/"."$_"."_2".".fastq";
       my $read="./temp_chr/"."$_".".fastq";
       system ("perl step2_sumCm.pl -ref $opt{ref} -1 $read1 -2 $read2 -read $read -project $project > $_.log 2> $_.log2");        
   }
}elsif(exists $opt{1} and exists $opt{2}){
   foreach (keys %chr){  
       my $project="./temp_chr/$_";
       my $read1="./temp_chr/"."$_"."_1".".fastq";
       my $read2="./temp_chr/"."$_"."_2".".fastq";
       system ("perl step2_sumCm.pl -ref $opt{ref} -1 $read1 -2 $read2 -project $project > $_.log 2> $_.log2");
   }     
}elsif(exists $opt{read}){
   foreach (keys %chr){
       my $project="./temp_chr/$_";
       my $read="./temp_chr/"."$_".".fastq"; 
       system ("perl step2_sumCm.pl -ref $opt{ref} -read $read -project $project > $_.log 2> $_.log2");
   }
}
system ("cat ./temp_chr/*.plus.status > ./temp_chr/$opt{project}.plus.status");
system ("cat ./temp_chr/*.minus.status > ./temp_chr/$opt{project}.minus.status");
print "Done!\n";


### Split BED file and return three hash, %readSE (unpaired), %readPE (paired).
### keys in hash are read name and values are mapped chromosome (SRR032638.1->chr8).
sub splitBED
{
my ($project,$suffix)=@_;
my %read;
my $file=$project.$suffix;
open IN, "$file" or die "$!";
while (<IN>){
   my @unit=split("\t",$_);
   if ($unit[3]=~/(.*)\/\d+/){
       $read{$1}=$unit[0];
   }else{
       $read{$unit[3]}=$unit[0];
   }
   $chr{$unit[0]}=1;
   my $outfile="./temp_chr/$unit[0]".$suffix;
   my $content=$_;
   writefile($outfile,$content);       
}
close IN;
return \%read;
}

### split fastq read 
sub splitread
{
my ($read,$pair,$readT,$readA)=@_; ### $read is a reference of @read. $pair==1 or 2 if paired, $pair==0 if unpaired
my $suffix;
if ($pair==0){
   $suffix=".fastq";
}else{
   $suffix="_".$pair.".fastq";
}
for(my $i=0;$i<@$read;$i++){
         print "$i\t$$read[$i]\n";
         open IN, "$$read[$i]" or die "$!";
             while (<IN>){
                  if ($_=~/^@(.*)$/){
                       my @unit=split(" ",$1);
                       my $head=$unit[0];
                       my $seq=<IN>;
                       my $head1=<IN>;
                       my $qual=<IN>;
                       my $content=$_.$seq.$head1.$qual;
                       if (exists $readT->{$head}){
                          my $file="./temp_chr/$readT->{$head}".$suffix;
                          writefile($file,$content);
                       }
                       if(exists $readA->{$head}){
                          my $file="./temp_chr/$readA->{$head}".$suffix;
                          writefile($file,$content);
                       }
                  }     
             }
         close IN;
}
}

### write fastq
sub writefile
{
my ($file,$content)=@_;
open OUT, ">>$file" or die "$!";
    print OUT "$content";
close OUT;
}


