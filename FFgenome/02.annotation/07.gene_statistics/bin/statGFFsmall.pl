#!/usr/bin/perl

use Getopt::Long;
my $help=<<USAGE;

perl statGFFgene.pl -gff gene.gff -fasta genome.fa > output &

USAGE


my %opt;
GetOptions(\%opt,"gff:s","fasta:s","help");

if ($opt{help} or keys %opt < 1 ){
   print "$help";
   exit();
}

my $genenum;   ### count total gene number=number of mRNA
my @genelen;   ### store gene length in array, including 5' UTR, 3' UTR, CDS, intron=mRNA end - mRNA start +1  
my @cdslen;    ### coding sequence length
my @exonlen;   ### coding sequence length=add CDS length
my @exonnum;   ### exon number in each gene
my @intronlen; ### every intron length
my @genegc;    ### gc content of gene, exon, intron, UTR
my @cdsgc;     ### gc of cds
my @exongc;    ### gc of exon
my @introngc;  ### gc of intron
my $single;    ### single exon genes

my $refgff=parseGFF($opt{gff});
my $refseq=getfastaseq($opt{fasta});
foreach (keys %$refgff){
     my $index=$_;
     $genenum++;
     my @unit=split("\n",$refgff->{$_});  
     my $sequence;  
     my $exonnum;
     my @exon;
     my $CDSlen;
     my $CDS;
     my $strand;
     foreach (@unit){
          my @array=split("\t",$_);
          if ($array[2] eq "mRNA"){
                #print "$index\n";
                $sequence=$array[0];
                $strand=$array[6];
                my $len=$array[4]-$array[3]+1;
                my $geneseq=substr($refseq->{$sequence},$array[3]-1,$len);
                my $genegc =estimateGC($geneseq);
                push (@genelen,$len);
                push (@genegc, $genegc);
          }elsif($array[2] eq "CDS"){
                $exonnum++;
                push (@exon,[@array]);
          }
     }
          @exon=sort {$a->[3] <=> $b->[3]} @exon;
          push (@exonnum,$exonnum);
          if ($exonnum > 1){
               for(my $i=0;$i<@exon-1;$i++){
                   #print "$exon[$i+1][3]\t$exon[$i][4]\t";
                   my $len=$exon[$i+1][3]-$exon[$i][4]+1;
                   #print "$len\n";
                   my $intronseq=substr($refseq->{$sequence},$exon[$i][4],$len);
                   my $introngc =estimateGC($intronseq);
                   push (@intronlen,$len);
                   push (@introngc, $introngc);
               }
          }
          if ($exonnum == 1){
               $single++;
          }
               for(my $i=0;$i<@exon;$i++){
                   my $len=$exon[$i][4]-$exon[$i][3]+1;
                   $CDSlen+=$len;
                   my $exonseq=substr($refseq->{$sequence},$exon[$i][3]-1,$len);
                   my $exongc =estimateGC($exonseq);
                   $CDS.=$exonseq;
                   push (@exonlen,$len);
                   push (@exongc,$exongc);
               }
               push (@cdslen,$CDSlen);
               if ($strand eq "-"){
                 $CDS=revcom($CDS);
               }
               my $CDSgc=estimateGC($CDS);
               push (@cdsgc,$CDSgc);
               #print "$CDSlen\n$CDS\n$CDSgc\n";
          
          
}

my ($meanglen,$medianglen)=statnum(@genelen);
my ($meanggc, $medianggc) =statnum(@genegc);
my ($meancgc, $mediancgc) =statnum(@cdsgc);
my ($meanclen,$medianclen)=statnum(@cdslen);
my ($meanenum,$medianenum)=statnum(@exonnum);
my ($meanelen,$medianelen)=statnum(@exonlen);
my ($meanegc, $medianegc) =statnum(@exongc);
my ($meanilen,$medianilen)=statnum(@intronlen);
my ($meanigc, $medianigc) =statnum(@introngc);

print "Total number of Gene:\t$genenum\n";
print "Gene Length (Mean/Median):\t$meanglen/$medianglen\n";
print "Gene GC (Mean/Median):\t$meanggc/$medianggc\n";
print "Coding sequence Length (Mean/Median):\t$meanclen/$medianclen\n";
print "Coding sequence GC (Mean/Median):\t$meancgc/$mediancgc\n";
print "Exon Number (Mean/Median):\t$meanenum/$medianenum\n";
print "Exon Length (Mean/Median):\t$meanelen/$medianelen\n";
print "Exon GC (Mean/Median):\t$meanegc/$medianegc\n";
print "Intron Length (Mean/Median):\t$meanilen/$medianilen\n";
print "Intron GC (Mean/Median):\t$meanigc/$medianigc\n";
print "Single Exon Gene:\t$single\n";

#########################sub functions###################################
sub revcom
{
my ($seq)=@_;
my $rev=reverse $seq;
$rev=~tr/ATGCatgc/TACGtacg/;
return $rev;
}

sub statnum
{
my (@num)=@_;
my ($total,$number,$min,$max,$mean,$median);
@num=sort {$a<=>$b} @num;
my $loop=0;
foreach  (@num) {
        $total+=$_;
        $loop++;
}
$number=@num;
$min=$num[0];
$max=$num[$number-1];
$mean=$total/$number;
$median=$num[int $number/2];
if ($mean > 10){
   $mean =sprintf("%d",$mean);
   $median=sprintf("%d",$median);
}elsif($mean > 1){
   $mean =sprintf("%.1f",$mean);
   $median =sprintf("%.1f",$median);
}else{
   $mean =sprintf("%.2f",$mean);
   $median =sprintf("%.2f",$median);
}
return ($mean,$median);
}


sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $seq;   ##Scaffold
my $id;    ##ID for element
my $record;##all line for this record
my $index; ##key, Seq_id
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/mRNA/){
        $seq=$unit[0];
        if ($unit[8]=~/ID=(.*);/ or $unit[8] =~/ID=(.*)/){
            $id=$1;
        }
        $record="$_\n";
        $index="$seq"."_"."$id";
        $hash{$index}=$record;
    }elsif($unit[0] eq $seq and $unit[8] =~ /Parent=$id/){
        $hash{$index}.="$_\n";
    }

}
close IN;

return \%hash;
}


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

sub estimateGC{
my ($seq)=@_;
my $length=length $seq;
#my $a=$seq=~tr/Aa/Aa/g;
my $n=$seq=~tr/Nn/Nn/;
my $len=$length-$n;
my $c=$seq=~tr/Cc/Cc/;
my $g=$seq=~tr/Gg/Gc/;
if ($len > 0){
  my $gc=($g+$c)/$len;
  return $gc;
}else{
  my $gc="Na";
  return $gc;
}
}


