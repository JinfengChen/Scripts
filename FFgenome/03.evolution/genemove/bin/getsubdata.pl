#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"genegff:s","tegff:s","fasta:s","refgene:s","qrygene:s","flanking:s","draw","project:s","help");


my $help=<<USAGE;
Can get subfregment of gene gff and te gff for one gene or two genes. And draw dotplot of two regions with annotation.
perl $0 --genegff --tegff --fasta --refgene --qrygene --flanking --draw --project Brd

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

$opt{project} ||="Brd";
$opt{flanking} ||=25000;

my $refgff=parseGFF($opt{genegff});
my $refchr=$refgff->{$opt{refgene}}->[0];
my $refstart=$refgff->{$opt{refgene}}->[1]-$opt{flanking};
my $refend  =$refgff->{$opt{refgene}}->[2]+$opt{flanking};
my $refhead="$opt{project}\_$refchr\_$refstart\_$refend";
print "Refhead:$refhead\n";
& getsubfasta($opt{fasta},$refhead,"+");
& getsubgff3($opt{genegff},$refhead,"+");
& getsubrepeat($opt{tegff},$refhead,"+");

my ($qrychr,$qrystart,$qryend,$qryhead);
if (defined $opt{qrygene}){
$qrychr=$refgff->{$opt{qrygene}}->[0];
$qrystart=$refgff->{$opt{qrygene}}->[1]-$opt{flanking};
$qryend  =$refgff->{$opt{qrygene}}->[2]+$opt{flanking};
$qryhead="$opt{project}\_$qrychr\_$qrystart\_$qryend";
print "Qryhead:$qryhead\n";
& getsubfasta($opt{fasta},$qryhead,"+");
& getsubgff3($opt{genegff},$qryhead,"+");
& getsubrepeat($opt{tegff},$qryhead,"+");
}

if ($opt{draw}){
`perl dotplotRegion.pl --refseq $refhead.fasta --qryseq $qryhead.fasta --refgenegff $refhead.gene.gff --qrygenegff $qryhead.gene.gff --project test`;
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
        $hash{$id}=[$unit[0],$unit[3],$unit[4],$unit[6]];
    }

}
close IN;
return \%hash;
}

#########
sub getsubrepeat
{
my ($gff,$head,$strand)=@_;
my @temp=split("_",$head);
my $chr=$temp[1];
my $start=$temp[2];
my $end  =$temp[3];
my $len=$end-$start+1;
#print "$chr\t$start\t$end\n";
if ($strand eq "+"){
   open IN1, "$gff" or die "$!";
   open OUT1, ">>$head.te.gff" or die "$!";
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
   open OUT1, ">>$head.te.gff" or die "$!";
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
   open OUT1, ">>$head.gene.gff" or die "$!";
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
   open OUT1, ">>$head.gene.gff" or die "$!";
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
close IN;
$/="\n";
return \%hash;
}
 
