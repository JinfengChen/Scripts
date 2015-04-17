#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"help");


my $help=<<USAGE;
This script is design to select BAC clone for FISH experiment. We intend to select these BAC clones, which have pair of BES hit on scaffold and have a repeat content less than 10%. Also, only single hit BES pairs are prefered.
Usage: perl BAC4FISH.pl Scaffold000001 > Scaffold000001.BAC4FISH
The result will be a list a clone, following with hit position on scaffold, repeat content. 
USAGE

if ($opt{help} or @ARGV < 1){
   print "$help\n";
   exit();
}

my %bes;
open IN, "../input/bes.length" or die "$!";
while (<IN>){
     chomp $_;
     my @unit=split("\t",$_);
     if ($unit[0]=~/OB__B(\w+)\.(\w+)/){
          if (exists $bes{$1}){
              $bes{$1}.="\t$2\t$unit[1]\t$unit[2]";
          }else{
              $bes{$1}="$2\t$unit[1]\t$unit[2]";
          }
     }
}
close IN;


my %ref;
my %qry;
my %ctg;
my %beshit;
my @bes2scaf;
my @scafbes;
my $flag;
open IN, "../input/FFversion2_BSS_results_super-scaffold.brachyanthaBES.filter.bss" or die "$!";
while (<IN>){
    chomp $_;
    my @unit=split(" ",$_);
    if ($_=~/#Hit table/){
       <IN>;
       $flag=1;
    }elsif ($flag){
      my $record=join("\t",@unit);
      if (exists $beshit{$unit[3]}){
        my $rev=substr($unit[0],8,1);
        #print "$rev\n";
        if ($rev eq "f"){
           my $f=$beshit{$unit[3]}->[0];
           $f++;
           $beshit{$unit[3]}->[0]=$f;
        }else{
           my $r=$beshit{$unit[3]}->[1];
           $r++;
           $beshit{$unit[3]}->[1]=$r;
        }
      }else{
        my $rev=substr($unit[0],8,1);
        if ($rev eq "f"){
           $beshit{$unit[3]}=[1,0];
        }else{
           $beshit{$unit[3]}=[0,1];
        }
      }
      if ($unit[5] eq $ARGV[0]){
          push (@bes2scaf,$record);
          push (@scafbes, $unit[3]);
      }
     }
}
close IN;


## store BES pair of specific scaffold that have only one hit on all scaffold
my %record2scaf;
foreach(@bes2scaf){
    my @unit=split("\t",$_);
    if ($beshit{$unit[3]}->[0] == 1 and $beshit{$unit[3]}->[1] == 1){
       $record2scaf{$unit[0]}=$_;
    }

}
my $counter=0;
my $sequence=Scaffold($ARGV[0]);  ## get the scaffold sequence
@scafbes=uniqarray(@scafbes);     ## kick duplicate clone name
print "Clone\tstart\tend\tlen\trepeat\tgap\ttotalN\n";
foreach(@scafbes){  
   my $forward=$_."f";
   my $reverse=$_."r";
   if (exists $record2scaf{$forward} and exists $record2scaf{$reverse}){
         $counter++;
         my @f=split(" ",$record2scaf{$forward});
         my @r=split(" ",$record2scaf{$reverse});
         my ($start,$end,$len);
         if ($f[1] eq "n"){
              $start=$f[11];
              $end  =$r[12];
              $len  =$end-$start+1;
         }else{
              $start=$r[11];
              $end  =$f[12];
              $len  =$end-$start+1;
         }
         my $cloneseq=substr($sequence,$start,$len);
         #my ($repeat,$gap)=RepeatContent($cloneseq); 
         my $repeat=RepeatContentN($cloneseq); 
         my $totalN=$repeat;
         if ($totalN <= 1 and $len >= 50000){
            print "$_\t$start\t$end\t$len\t$repeat\t$gap\t$totalN\n";
         }
   }
}
print "Totally pair end:\t$counter\n";


## kick duplicated record in array
sub uniqarray
{
my (@array)=@_;
my %hash;
foreach (@array){
   $hash{$_}=1;
}
my @temp=keys %hash;
return @temp;
}


##Estimate repeat content for DNA sequence
##We treat small case of charactar as repeat, and N as gap. 
##Report both of them.
sub RepeatContent
{
my ($seq)=@_;
my $length=length $seq;
my $g=$seq=~tr/Nn/Nn/;
my $n=$seq=~tr/atcg/atcg/;
my $repeat=$n/$length;
my $gap=$g/$length;
return ($repeat,$gap);
}

sub RepeatContentN
{
my ($seq)=@_;
my $length=length $seq;
my $n=$seq=~tr/Nn/Nn/;
my $repeat=$n/$length;
return $repeat;
}



## get scaffold sequence from repeat masked file of super-scaffold.fa
sub Scaffold
{
$/=">";
my ($name)=@_;
my $seq;
open IN, "../input/FFversion2/super-scaffold.allmask.fa" or die "$!";
while (<IN>){
    chomp $_;
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $head=shift @unit;
    my @array=split(" ",$head);
    $head=$array[0];
    my $sequence=join("",@unit);
    if ($head eq $name){
        $seq=$sequence;
        last;
    }
}
close IN;
$/="\n";
return $seq;
}



