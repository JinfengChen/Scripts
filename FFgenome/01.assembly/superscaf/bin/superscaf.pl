#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"fasta:s","link:s","cut:s","project:s","help");


my $help=<<USAGE;
Link super-scaffold from fpc based manual check table.
The scaffold sequence that mapped to the fpc physical map is also called super-scaffold, which construted by BES scaffolding. 
After link, we create two new sequence file, one for super-scaffold and the other for chromosome.
-fasta:   scaffold sequence
-link:    link table
-cut:   contain split point table for scaffold which are wrong assembled, optional.
-project: project name used to name output sequence file

Example:
Chr	Super	Scaffold	Strand	Evidence
chr1	1	Scaffold000152	plus	rice/fpc
chr1	1	Scaffold000414	minus	fpc

USAGE

if ($opt{help} or keys %opt < 1){
   print "$help\n";
   exit ();
}

my $scafseq=getfastaseq($opt{fasta});
our $totallen;
##########split scaffold according the cut point, which are identified to be misassembly###############################
if (exists $opt{cut}){
   open IN, "$opt{cut}" or die "$!";
      while (<IN>){
         next if (length $_ < 2);
         my @unit=split("\t",$_);
         if (exists $scafseq->{$unit[0]}){
             my $aname=$unit[0]."a";
             my $bname=$unit[0]."b";
             my $aseq =substr($scafseq->{$unit[0]},0,$unit[1]);
             my $bseq =substr($scafseq->{$unit[0]},$unit[1]);
             my $newaseq=trim3N($aseq);
             my $newbseq=trim5N($bseq); 
             $scafseq->{$aname}=$newaseq;
             $scafseq->{$bname}=$newbseq;
             delete $scafseq->{$unit[0]};
         }
      }
   close IN;
}
################################################################################################################

###############read data from link table and join the scaffold################################################
my $n="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
my %chr;
my %super;
open IN, "$opt{link}" or die "$!";
<IN>;
while (<IN>){
    my @unit=split("\t",$_);
    if ($unit[3] eq "plus"){
       if (exists $scafseq->{$unit[2]}){ 
          if (exists $chr{$unit[0]}){
              $chr{$unit[0]}.=$n.$scafseq->{$unit[2]};
          }else{
              $chr{$unit[0]}=$scafseq->{$unit[2]};
          }
          if (exists $super{$unit[1]}){
              $super{$unit[1]}.=$n.$scafseq->{$unit[2]};
          }else{
              $super{$unit[1]}=$scafseq->{$unit[2]};
          }
          delete $scafseq->{$unit[2]};
       }else{
          print "Can not found sequence of $unit[2]\n";
       }
    }else{
       if (exists $scafseq->{$unit[2]}){
          my $recseq=revcom($scafseq->{$unit[2]});
          if (exists $chr{$unit[0]}){
              $chr{$unit[0]}.=$n.$recseq;
          }else{
              $chr{$unit[0]}=$recseq;
          }
          if (exists $super{$unit[1]}){
              $super{$unit[1]}.=$n.$recseq;
          }else{
              $super{$unit[1]}=$recseq;
          }
          delete $scafseq->{$unit[2]};
       }else{
          print "Can not found sequence of $unit[2]\n";
       }
    }

}
close IN;
################# write chromosome and super into files###########################################
open OUT, ">$opt{project}.chr.fa" or die "$!";
foreach (sort keys %chr){
    my $fseq=formatseq($chr{$_},50);
    print OUT ">$_\n$fseq\n";
}
close OUT;

open OUT, ">$opt{project}.super.fa" or die "$!";
foreach (sort {$a <=> $b} keys %super){
    my $fseq=formatseq($super{$_},50);
    printf OUT (">Scaffold%06d\n",$_);
    print OUT "$fseq\n";
}
close OUT;

my $leftlen;
my $leftscafn=36;
open OUT, ">$opt{project}.unmaped.fa" or die "$!";
foreach (sort keys %$scafseq){
   $leftscafn++;
   $leftlen+=length $scafseq->{$_};
   my $fseq=formatseq($scafseq->{$_},50);
   my $leftscaf=sprintf("Scaffold%06d",$leftscafn);
   print OUT ">$leftscaf\n$fseq\n";
}
close OUT;

my $rate=$leftlen/$totallen;
print "Left length: $leftlen\n";
print "Percent: $rate\n";
##############################################################################################################


=pod
foreach (keys %$scafseq){
    my $seq=formatseq($scafseq->{$_},100);
    print "$seq\n";
}
for (my $i=1;$i<=10;$i++){
    printf ("Scaffold%08d\n",$i);
    
}
=cut


sub trim5N
{
my ($seq)=@_;
my $revseq=reverse $seq;
my $trimseq=trim3N($revseq);
my $newseq=reverse $trimseq;
return $newseq;
}

sub trim3N
{
my ($seq)=@_;
my $do=1;
while($do==1){
  my $last=substr($seq,-1);
  if ($last =~/N/i){
     chop $seq
  }else{
     $do=0;
  }
}
return $seq;
}



sub formatseq
{
### format a single line sequence into lines with user specific length
my ($seq,$step)=@_;
my $length=length $seq;
my $run=int ($length/$step);
my $newseq;
for(my $i=0;$i<=$run;$i++){
   my $start=$i*$step;
   my $line=substr($seq,$start,$step);
   $newseq.="$line\n";
}
return $newseq;
}


sub revcom
{
my ($seq)=@_;
my $rev=reverse $seq;
$rev=~tr/ATGCatgc/TACGtacg/;
return $rev;
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
    $totallen+=length $seq;
    $hash{$head}=$seq;
}
$/="\n";
print "Total length: $totallen\n";
return \%hash;
}

