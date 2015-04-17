#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"fasta:s","gene:s","te:s","tandem:s","help");


my $help=<<USAGE;
This script is designed to get tandem repeat region sequence and gff annotation.

Run: perl getTandemRegion.pl -fasta genome.fa -gene gene.gff -te te.gff -tandem test.tandem_repeat_gene.txt
-fasta: fasta format of genome sequence 
-gene:  gff format gene annotation
-te:    gff format te annotation
-tandem:tandem_repeat_gene.txt, one line per cluster
Example tandem_repeat_gene.txt:
OBR_GLEAN_10020666 OBR_GLEAN_10020657
OBR_GLEAN_10027586 OBR_GLEAN_10027584
OBR_GLEAN_10007950 OBR_GLEAN_10007942 OBR_GLEAN_10007949



USAGE

if ($opt{help} or keys %opt < 1){
       print "$help\n";
       exit();
}


my ($refgff,$refstart,$refname)=parseGFF($opt{gene});
my $refseq=getfastaseq($opt{fasta});

my $cluster;
open IN, "$opt{tandem}" or die "$!";
open OUT, ">Ob_tandemRegion.fa" or die "$!";
open TAN, ">TandemInfTable.txt" or die "$!";
while (<IN>){
    chomp $_;
    next if (length $_ < 2);
    $cluster++;
    my @unit=split(" ",$_);
    my @startp;
    my %ref;
    my $seqname=$refname->{$unit[0]};
    my $refnum=1;
    my $member=@unit;
    foreach (@unit){
       push (@startp,$refstart->{$_});
       if ($seqname ne $refname->{$_}){
           $refnum++;
       }
    }
    unless ($refnum == 1){
        print "Multi reference for cluster\n";
        next;
    }
    @startp=sort {$a <=> $b} @startp;
    my $start=shift @startp;
    my $end  =pop @startp;
    my $len  =$end -$start +1;
    #print "Cluster$cluster: $seqname\t$start\t$end\t$len\n";
    my $seq=$refseq->{$seqname};
    my $seqlen=length $seq;
    if ($start-150000 <= 0){
        $start=0;
    }else{
        $start=$start-150000;
    }
    if ($end+150000 >= $seqlen-1){
        $end=$seqlen-1;
    }else{
        $end=$end+150000;
    }
    my $finallen=$end-$start+1;
    print "Cluster$cluster: $member\t$seqname\t$start\t$end\t$finallen\n";
    my $tandemseq=substr($seq,$start,$finallen);
    my $head="OB_".$seqname."_".$start."_".$end;
    if ($member >= 5){
        my $m=join("\t",@unit);
        print TAN "Cluster$cluster\t$member\t$head\t$m\n";
        print OUT ">$head\n$tandemseq\n";
        getsubGFF($opt{gene},$head);
        getsubGFF($opt{te},$head);
    }
    
}
close IN;
close OUT;
close TAN;
#####################sub function blocks#############################

sub getsubGFF
{
my ($gff,$head)=@_;
my @temp=split("_",$head);
my $chr=$temp[1];
my $start=$temp[2];
my $end  =$temp[3];
#print "$chr\t$start\t$end\n";
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
return 1;
}


sub parseGFF
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of id
my %start; ##hash to store gene start position
my %ref;
my $seq;   ##Scaffold
my $id;    ##ID for element
my $record;##all line for this record
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
        $hash{$id}=$record;
        $start{$id}=$unit[3];
        $ref{$id}=$seq;
    }elsif($unit[0] eq $seq and $unit[8] =~ /Parent=$id/){
        $hash{$id}.="$_\n";
    }

}
close IN;

return (\%hash,\%start,\%ref);
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

