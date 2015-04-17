#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"qryhead:s","qrygff3:s","qryrepeat:s","qryfasta:s","qRNAseq:s","qmC:s","qH3K4:s","qH3K9:s","rRNAseq:s","rmC:s","rH3K4:s","rH3K9:s","refhead:s","refgff3:s","refrepeat:s","reffasta:s","outdir:s","help");


my $help=<<USAGE;
Prepare ACT files and BED alignment file for a region.
These files are then used to draw comparative figure and signal graph using drawACTgraph.pl
--qryhead :  query name that used to pass position information,OB_chr01_16397572_16441992
--qrygff3 :  dir of gff3 file for each chromosome
--qryfasta:  dir of fasta file for each chromosome
--qryrepeat: dir of gff file of repeat
--qRNAseq :  dir of bed file of RNAseq
--qmC     :  dir of bed file of methylation
--qH3K4   :  dir of bed file of H3K4me2
--qH3K9   :  dir of bed file of H3K9me2
--rRNAseq :
--rmC     :
--rH3K4   :
--rH3K9   :
--refhead :
--refgff3 :
--reffasta:
--refrepeat:
--outdir  :  output dir
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

$skipTE=1;
$opt{outdir} ||=".";

& createACT($opt{qryhead},$opt{refhead},"$opt{outdir}/$opt{refhead}");
###get sub bed ##########
if (exists $opt{qRNAseq}){
    & getsubbed($opt{qryhead},$opt{qRNAseq},"$opt{outdir}/$opt{refhead}","RNAseq");
}
if(exists $opt{qmC}){  
    & getsubbed($opt{qryhead},$opt{qmC},"$opt{outdir}/$opt{refhead}","mC");
}
if(exists $opt{qH3K4}){
    & getsubbed($opt{qryhead},$opt{qH3K4},"$opt{outdir}/$opt{refhead}","H3K4");
}
if(exists $opt{qH3K9}){
    & getsubbed($opt{qryhead},$opt{qH3K9},"$opt{outdir}/$opt{refhead}","H3K9");
}
if (exists $opt{rRNAseq}){
    & getsubbed($opt{refhead},$opt{rRNAseq},"$opt{outdir}/$opt{refhead}","RNAseq");
}
if(exists $opt{rmC}){
    & getsubbed($opt{refhead},$opt{rmC},"$opt{outdir}/$opt{refhead}","mC");
}
if(exists $opt{rH3K4}){
    & getsubbed($opt{refhead},$opt{rH3K4},"$opt{outdir}/$opt{refhead}","H3K4");
}
if(exists $opt{rH3K9}){
    & getsubbed($opt{refhead},$opt{rH3K9},"$opt{outdir}/$opt{refhead}","H3K9");
}
#system ("/home/jfchen/FFproject/tools/BEDTools/bin/windowBed -a a.bed -b b.bed -w 0 -u")

############# sub functions ############ 
sub getsubbed
{
my ($head,$bed,$dir,$title)=@_;
my @inf=split("_",$head);
open TEMP, ">temp.bed" or die "$!";
print TEMP "$inf[1]\t$inf[2]\t$inf[3]\n";
close TEMP;
my $beda="$bed/$inf[1]";
my $bedb="temp.bed";
my $bedout="$head"."."."$title";
print "$bedout\n";
`/home/jfchen/FFproject/tools/BEDTools/bin/windowBed -a $beda -b $bedb -w 1 -u > $bedout`;
open BED1, "$bedout" or die "$!";
open BED2, ">$bedout.bed" or die "$!";
while(<BED1>){
   chomp $_;
   next if ($_ eq "");
   my @array=split("\t",$_);
   if ($array[1] >= $inf[2]){
       $array[1]=$array[1]-$inf[2];
       $array[2]=$array[2]-$inf[2];
       my $temp=join("\t",@array);
       print BED2 "$temp\n";
   }
}
close BED2;
close BED1;
`mv $bedout* $dir/`;
}


#######

sub createACT
{
my ($qryhead,$refhead,$dir)=@_;
`mkdir $dir` unless (-f $dir);
my @qryinf=split("_",$qryhead);
my @refinf=split("_",$refhead);
getsubgff3("$opt{qrygff3}/$qryinf[1]",$qryhead,"+");
getsubgff3("$opt{refgff3}/$refinf[1]",$refhead,"+");
getsubfasta("$opt{qryfasta}/$qryinf[1]",$qryhead,"+");
getsubfasta("$opt{reffasta}/$refinf[1]",$refhead,"+");
getsubrepeat("$opt{qryrepeat}/$qryinf[1]",,$qryhead,"+") if (-f "$opt{qryrepeat}/$qryinf[1]");
getsubrepeat("$opt{refrepeat}/$refinf[1]",,$refhead,"+") if (-f "$opt{refrepeat}/$refinf[1]");
`perl GFF2embl.pl -gff $qryhead.gff -embl $qryhead.embl -fasta $qryhead.fasta`;
`perl gffrepeat2embl.pl -repeat $qryhead.repeat -embl $qryhead.embl -title $qryhead` if (-f "$opt{qryrepeat}/$qryinf[1]");
`mv $qryhead.embl $qryhead.gene` if (-f "$opt{qryrepeat}/$qryinf[1]");
`mv $qryhead.merge $qryhead.embl` if (-f "$opt{qryrepeat}/$qryinf[1]");
`perl GFF2embl.pl -gff $refhead.gff -embl $refhead.embl -fasta $refhead.fasta`;
`perl gffrepeat2embl.pl -repeat $refhead.repeat -embl $refhead.embl -title $refhead` if (-f "$opt{refrepeat}/$refinf[1]");
`mv $refhead.embl $refhead.gene` if (-f "$opt{refrepeat}/$refinf[1]");
`mv $refhead.merge $refhead.embl` if (-f "$opt{refrepeat}/$refinf[1]");
#`perl runblast2seq.pl`;
#`perl run2act.pl`;
#`rm *.fasta.n* *.blast *.temp *.gff *.fasta *.gene *.repeat`;
#`rm *.fasta.n* *.blast *.temp *.gff *.gene *.repeat`;
`mv $qryhead* $dir/`;
`mv $refhead* $dir/`;
}
#####################
sub getsubrepeat
{
my ($gff,$head,$strand)=@_;
my @temp=split("_",$head);
my $chr=$temp[1];
my $start=$temp[2] > $temp[3] ? $temp[3] : $temp[2];
my $end  =$temp[2] > $temp[3] ? $temp[2] : $temp[3];
my $len=$end-$start+1;
#print "$chr\t$start\t$end\n";
if ($strand eq "+"){
   open IN1, "$gff" or die "$!";
   open OUT1, ">>$head.repeat" or die "$!";
   while (<IN1>){
     my @unit=split("\t",$_);
     if ($unit[0] eq $chr and $unit[3] >= $start and $unit[4] <= $end){
       $unit[3]=$unit[3]-$start;
       $unit[4]=$unit[4]-$start;
       my $line=join("\t",@unit);
       print OUT1 "$line";
     }

   }
   close OUT1;
   close IN1;
}else{
   open IN1, "$gff" or die "$!";
   open OUT1, ">>$head.repeat" or die "$!";
   while (<IN1>){
     my @unit=split("\t",$_);
     if ($unit[0] eq $chr and $unit[3] >= $start and $unit[4] <= $end){
        my $tempend   =$len-($unit[3]-$start);
        my $tempstart =$len-($unit[4]-$start);
        $unit[3]=$tempstart;
        $unit[4]=$tempend;
        if ($unit[6] eq "+"){
            $unit[6] = "-";
        }else{
            $unit[6] = "+";
        }
        my $line=join("\t",@unit);
        print OUT1 "$line";
     }
   }
   close OUT1;
   close IN1;
}
}

#####################
sub getsubgff3
{
my ($gff,$head,$strand)=@_;
my @temp=split("_",$head);
my $chr=$temp[1];
my $start=$temp[2] > $temp[3] ? $temp[3] : $temp[2];
my $end  =$temp[2] > $temp[3] ? $temp[2] : $temp[3];
my $len=$end-$start+1;
print "$gff\t$chr\t$start\t$end\n";
if ($strand eq "+"){
   open IN1, "$gff" or die "$!";
   open OUT1, ">>$head.gff" or die "$!";
   while (<IN1>){
     my @unit=split("\t",$_);
     if ($unit[0] eq $chr and $unit[3] >= $start and $unit[4] <= $end){ 
       $unit[3]=$unit[3]-$start;
       $unit[4]=$unit[4]-$start;
       my $line=join("\t",@unit);
       print OUT1 "$line";
     } 

   }
   close OUT1;
   close IN1;
}else{
   open IN1, "$gff" or die "$!";
   open OUT1, ">>$head.gff" or die "$!";
   while (<IN1>){
     my @unit=split("\t",$_);
     if ($unit[0] eq $chr and $unit[3] >= $start and $unit[4] <= $end){ 
        my $tempend   =$len-($unit[3]-$start);
        my $tempstart =$len-($unit[4]-$start); 
        $unit[3]=$tempstart;
        $unit[4]=$tempend;
        if ($unit[6] eq "+"){
            $unit[6] = "-";
        }else{
            $unit[6] = "+";
        }
        my $line=join("\t",@unit);
        print OUT1 "$line";
     }
   }
   close OUT1;
   close IN1;
}
}
####################
sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;
my $id;    ##ID for element
my $record;##all line for this record
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*?);/){
            $id=$1;
        }
        $record="$_\n";
        $hash{$id}=$record;
    }elsif($unit[0] eq $seq and $unit[8] =~ /Parent=$id/){
        $hash{$id}.="$_\n";
    }

}
close IN;
return \%hash;
}

################
sub getsubfasta
{
my ($fasta,$head,$strand)=@_;
my @temp=split("_",$head);
my $chr=$temp[1];
my $start=$temp[2] > $temp[3] ? $temp[3] : $temp[2];
my $end  =$temp[2] > $temp[3] ? $temp[2] : $temp[3];
my $len=$end-$start+1;
my $refseq=getfastaseq($fasta);
if ($strand eq "+"){
   if (exists $refseq->{$chr}){
       my $subseq=substr($refseq->{$chr},$start,$len);
       open FAS, ">$head.fasta" or die "$!"; 
             print FAS ">$head\n$subseq\n";
       close FAS;  
   }else{
       print "$chr can not found in $fasta\n";
   }
}else{
   if (exists $refseq->{$chr}){
       my $subseq=substr($refseq->{$chr},$start,$len);
       $subseqrec=revcom($subseq);
       open FAS, ">$head.fasta" or die "$!";
             print FAS ">$head rec\n$subseqrec\n";
       close FAS;
   }else{
       print "$chr can not found in $fasta\n";
   }   
}
}


################
sub getfastaseq
{
$/=">";
my %hash;
my ($file)=@_;
open IN,"$file" or die "$!";
while (<IN>){
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $temp=shift @unit;
    my @temp1=split(" ",$temp);
    my $head=$temp1[0];
    my $seq=join("\n",@unit);
    $seq=~s/\>//g;
    $seq=~s/\n//g;
    $seq=~s/\s//g;
    #print "$head\n";
    $hash{$head}=$seq;
}
$/="\n";
return \%hash;
}


