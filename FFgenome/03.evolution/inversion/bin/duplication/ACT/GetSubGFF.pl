#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"qryhead:s","qrygff3:s","qryfasta:s","project:s","help");


my $help=<<USAGE;
Get sub gff and fasta file according to qryhead position
--qryhead :  query name that used to pass position information,OB_chr01_16397572_16441992
--qrygff3 :  gff3 file for each chromosome
--qryfasta:
--project:
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

& createACT($opt{qryhead});

sub createACT
{
my ($qryhead)=@_;
my @qryinf=split("_",$qryhead);
getsubgff3("$opt{qrygff3}",$qryhead,"+");
getsubfasta("$opt{qryfasta}",$qryhead,"+");
`perl GFF2embl.pl -gff $qryhead.gff -embl $qryhead.embl -fasta $qryhead.fasta`;
}

sub getsubgff3
{
my ($gff,$head,$strand)=@_;
my @temp=split("_",$head);
my $chr=$temp[1];
my $start=$temp[2];
my $end  =$temp[3];
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
my $start=$temp[2];
my $end  =$temp[3];
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


