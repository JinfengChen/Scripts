#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"--maf:s","verbose","help");


my $help=<<USAGE;
perl $0 --maf
Summary align and gap length in maf alignment.
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $alignlen=0;
my $gaplen=0; # alignlen should be gap free alignmeng len
my @array;
open IN, "$opt{maf}" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    next if ($_=~/^#/);
    my $line=$_;
    if ($line=~/^a/){ # every aligned segment
       my $refseq=<IN>;
       my $ref=mafline($refseq);
       my $qryseq=<IN>;
       my $qry=mafline($qryseq);
       my ($align,$gap)=align($ref->[2],$qry->[2]);
       $alignlen+=$align;
       $gaplen+=$gap;
    }
}
close IN;

print "ALIGN: $alignlen\n";
print "GAP: $gaplen\n";

####################################
sub align
{
my ($copy,$check)=@_;
my $align=0;
my $gap=0;
my @ref=split("",$copy);
my @qry=split("",$check);
for (my $i=0;$i<@ref;$i++){
   if ($ref[$i]=~/[atcg]/i){
      $ref[$i]=~tr/atcg/ATCG/;
      $qry[$i]=~tr/atcg/ATCG/;
      if ($qry[$i] eq $ref[$i]){
         $align++;
      }else{
         $gap++;
         print "$i\t$qry[$i]\t$ref[$i]\n" if $opt{verbose};
      }
   }
}
return ($align,$gap);
}

###
sub mafline
{
my ($line)=@_;
my @unit=split(" ",$line);
my $start=$unit[2];
my $end  =$unit[2]+$unit[3]-1;
my $seq  =$unit[6]; #have gap
my @inf=($start,$end,$seq);
return \@inf;
}


