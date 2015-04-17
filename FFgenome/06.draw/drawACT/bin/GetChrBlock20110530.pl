#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"align:s","qrygff:s","refgff:s","help");


my $help=<<USAGE;
--align:   dir of chr align from mcscan/DAG
--qrygff:  gff for mcscan
--refgff:  gff for mcscan
Run: perl $0 --align ../data/chralign/ --qrygff ../data/Ob.gff --refgff ../data/Os.gff > blockinf
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my $qrygff=gff($opt{qrygff});
my $refgff=gff($opt{refgff}); 

my @align=glob("$opt{align}/chr*.align");
foreach(@align){
   my $chr;
   if ($_=~/\/(chr\d+)\.align/){
      $chr=$1;
   }
   print "$chr\n";
   my $qrychrgff=$qrygff->{$chr};
   my $refchrgff=$refgff->{$chr}; 
   my $refalign=parseDAG($_);
   ###sort block based on the position of first gene in rice genome
   my %blocksort;
   my %blockinf;
   foreach(keys %$refalign){
       my $blockn=$_;
       my $block=$refalign->{$_};
       shift @$block;
       my $startgene=$block->[0]->[1];
       my $endgene=$block->[@$block-1]->[1];
       my $blockstart=$refchrgff->{$startgene}->[0];
       my $blockend=$refchrgff->{$endgene}->[0];
       $blocksort{$_}=$blockstart;
       $blockinf{$_}=[$startgene,$blockstart,$endgene,$blockend];
       #print "$blockn\t$startgene\t$blockstart\t$endgene\t$blockend\n";
   }
   my @sortblock=sort {$blocksort{$a} <=> $blocksort{$b}} keys %blocksort;
   my $startref=0;
   my $startqry=0;
   foreach(@sortblock){
       #print "$_\t$blockinf{$_}->[0]\t$blockinf{$_}->[1]\t$blockinf{$_}->[2]\t$blockinf{$_}->[3]\n";
       my $block=$refalign->{$_};
=pod
       unless (defined $startref and defined $startqry){
           my $startrefgene=$block->[0]->[1];
           my $startqrygene=$block->[0]->[0];
           $startref=$refchrgff->{$startrefgene}->[0];
           $startqry=$qrychrgff->{$startqrygene}->[0];
       }
=cut
############# start block and print out head if there was a block interal
       my $temprefgene=$block->[0]->[1];
       my $tempqrygene=$block->[0]->[0];
       $temprefstart=$refchrgff->{$temprefgene}->[0];
       $tempqrystart=$qrychrgff->{$tempqrygene}->[0];
       $temprefend  =$refchrgff->{$temprefgene}->[1]+1000;
       $tempqryend  =$qrychrgff->{$tempqrygene}->[1]+1000;
       if ($temprefstart-$startref > 10000){
           my $refhead="OS"."_".$chr."_".$startref."_".$temprefend;
           my $qryhead="OB"."_".$chr."_".$startqry."_".$tempqryend;
           print "$qryhead\t$refhead\n";
           $startref=$temprefstart-1000 >= 0 ? $temprefstart-1000:0;
           $startqry=$tempqrystart-1000 >= 0 ? $tempqrystart-1000:0;
       }else{
           $startref=$temprefstart-1000 >= 0 ? $temprefstart-1000:0;
           $startqry=$tempqrystart-1000 >= 0 ? $tempqrystart-1000:0;
       }
######################################################################
       my $count;
       foreach(@$block){
           $count++;
           my $qrygene=$_->[0];
           my $refgene=$_->[1];
           my $qrystart=$qrychrgff->{$qrygene}->[0];
           my $refstart=$refchrgff->{$refgene}->[0];
           my $qryend  =$qrychrgff->{$qrygene}->[1]+1000;
           my $refend  =$refchrgff->{$refgene}->[1]+1000;
           if ($refstart-$startref >= 4000000){
               my $refhead="OS"."_".$chr."_".$startref."_".$refend;
               my $qryhead="OB"."_".$chr."_".$startqry."_".$qryend;  
               print "$qryhead\t$refhead\n";
               $startref=$refstart-1000;
               $startqry=$qrystart-1000;           
           }
           if ($count == @$block){

               my $refhead="OS"."_".$chr."_".$startref."_".$refend;
               my $qryhead="OB"."_".$chr."_".$startqry."_".$qryend;
               print "$qryhead\t$refhead\n";
               $startref=$refstart-1000;
               $startqry=$qrystart-1000;
           }
       }
   }
}

sub parseDAG
{
my ($align)=@_;
my %hash;
my $block;
my $qry;
my $ref;
open IN, "$align" or die "$!";
while(<IN>){
    chomp $_;
    if ($_=~/\#\# alignment\s+(\w+)\s+vs\.\s+(\w+)\s+Alignment\s+\#(\d+)\s+score/) {
    #if ($_=~/\#\# Alignment (\d+)\: score=.* e_value=.* N=\d+ (\w+)\&(\w+) (\w+)/) {
        $block=$3;
        $qry=$1;
        $ref=$2;
        #print "$block\t$qry\t$ref\n";
    }elsif($_=~/\#\# alignment\s+(\w+)\s+vs\.\s+(\w+)\s+.*?Alignment\s+\#(\d+)\s+score/){
        $block="R".$3;
        $qry=$1;
        $ref=$2;
    }elsif($_=~/$ref/){
        my @unit=split("\t",$_);
        #print "1:$unit[1]\t2:$unit[2]\n";
        if (exists $hash{$block}) {
           my $reftemp=$hash{$block};
           push (@$reftemp,[$unit[1],$unit[5]]);
           $hash{$block}=$reftemp;
        }else{
           my @temp;
           my @unit=split("\t",$_);
           push (@temp,[$qry,$ref]);
           push (@temp,[$unit[1],$unit[5]]);
           $hash{$block}=\@temp;
        }
     }
}
close IN;
return \%hash;
}


####################################################
sub parsemcscan
{
my ($align)=@_;
my %hash;
my $block;
my $qry;
my $ref;
open IN, "$align" or die "$!";
while(<IN>){
    chomp $_;
    if ($_=~/\#\# Alignment (\d+)\: score=.* e_value=.* N=\d+ (\w+)\&(\w+) (\w+)/) {
        $block=$1;
        $qry=$2;
        $ref=$3;
        #print "$block\t$qry\t$ref\n";
    }elsif($_=~/\s*$block\-\s*\d+\:/){
        my @unit=split("\t",$_);
        #print "1:$unit[1]\t2:$unit[2]\n";
        if (exists $hash{$block}) {
           my $reftemp=$hash{$block};
           push (@$reftemp,[$unit[1],$unit[2]]);
           $hash{$block}=$reftemp;
        }else{
           my @temp;
           my @unit=split("\t",$_);
           push (@temp,[$qry,$ref]);
           push (@temp,[$unit[1],$unit[2]]);
           $hash{$block}=\@temp;
        }
     }
}
close IN;
return \%hash;
}


sub gff
{
###deal with gff for mcscan
###Os01	Os01t001010	1023020	1200220
my ($gff)=@_;
my %hash;
open IN, "$gff" or die "$!";
while(<IN>){
        chomp $_;
        next if ($_ eq "");
        my @unit=split("\t",$_);
        if ($unit[0]=~/Ob/){
           $unit[0]=~s/Ob/chr/;
        }elsif($unit[0]=~/Os/){
           $unit[0]=~s/Os/chr/;
        }
        #print "$unit[0]\n";
        if (exists $hash{$unit[0]}) {
        my $reftemp=$hash{$unit[0]};
            $reftemp->{$unit[1]}=[$unit[2],$unit[3]];
            $hash{$unit[0]}=$reftemp;
        }else{
            my %temp;
            $temp{$unit[1]}=[$unit[2],$unit[3]];
            $hash{$unit[0]}=\%temp;
        }
}
close IN;
return \%hash;
}

