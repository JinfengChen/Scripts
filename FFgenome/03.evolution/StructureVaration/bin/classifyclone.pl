#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"blasttable:s","inf:s","project:s","help");


my $help=<<USAGE;
perl $0 --blasttable --inf
--inf FPC infromation of BES and clone.
--blasttable blasttable
--project: 
USAGE

if ($opt{help} or keys %opt < 1){
   print "$help\n";
   exit ();
}

my $fpcinf=readfpcinf($opt{inf});
my $title=$1 if ($opt{blasttable}=~/(.*)\.blasttable/);
$opt{project}||=$title;

my %rank=(
    1 => "TWO_OF_TWO_HIT_OPPDIR_TRUSTED",
    2 => "TWO_OF_TWO_HIT_OPPDIR_SHORT",
    3 => "TWO_OF_TWO_HIT_OPPDIR_LONG",
    4 => "TWO_OF_TWO_HIT_SAMEDIR_TRUSTED",
    5 => "TWO_OF_TWO_HIT_SAMEDIR_SHORT",
    6 => "TWO_OF_TWO_HIT_SAMEDIR_LONG",
    7 => "TWO_OF_TWO_HIT_OPPDIR_TOO_LONG",
    8 => "TWO_OF_TWO_HIT_SAMEDIR_TOO_LONG",
    9 => "TWO_OF_TWO_HIT_DIFF_CHROM",
    10=> "ONE_OF_TWO_READ_HIT"
);

my %valid=(
    1 => 1,
    2 => 0,
    3 => 0,
    4 => 0,
    5 => 0,
    6 => 0,
    7 => 0,
    8 => 0,
    9 => -1,
    10=> -1
);

my $min=70000;
my $max=200000;
my $maxlimit=500000;
####
my %summary;
my $clone=statblast("$opt{blasttable}");
open OUT, ">$title.classify" or die "$!";
foreach my $c (sort keys %$clone){
   #print "$bes\t$stat->{$bes}\n";
   my $class=classify($clone->{$c},$min,$max);
   my $crank=$rank{$class->[0]};
   $summary{$class->[0]}++;
   ##clone	rank	rankname	valid/discordant\n
   print OUT ">$c\n$class->[0]\t$crank\t$valid{$class->[0]}\n$class->[1]\n";
}
close OUT;
##print summary 
open OUT, ">$title.rank.summary" or die "$!";
foreach my $r (sort {$a <=> $b} keys %rank){
   $rcount=$summary{$r} > 0 ? $summary{$r} : 0;
   print OUT "$r\t$rank{$r}\t$rcount\n";
}
close OUT;



########################################
#Query_id        Query_length    Query_start     Query_end       Subject_id      Subject_length  Subject_start   Subject_end     Identity        Positive        Gap     Align_length    Score   E_value Query_annotation        Subject_annotation
sub classify
{
my ($clonehit,$min,$max,$maxlimit)=@_;
my $return;
if (!exists $clonehit->{0} or !exists $clonehit->{1}){
   my @hit= exists $clonehit->{0} > 0 ? @{$clonehit->{0}} : @{$clonehit->{1}};
   @hit=sort {$b->[12] <=> $a->[12]} @hit;
   #print "$hit[0]->[0]\t$hit[0]->[1]\n";
   my $line=join("\t",@{$hit[0]});
   $return=[10,$line];
}else{
   my $tmp=$clonehit->{1}->[0]->[0];
   #print "Temp:$tmp\n";
   my @hit; ### record all pair rank and hit information
   for(my $i=0;$i<@{$clonehit->{1}};$i++){
      for(my $j=0;$j<@{$clonehit->{0}};$j++){
         ## bes1
         my $chr1=$clonehit->{1}->[$i]->[4];
         my $strand1=$clonehit->{1}->[$i]->[7]-$clonehit->{1}->[$i]->[6] > 0 ? 1 : 0;
         my $start1 =$strand1 > 0 ? $clonehit->{1}->[$i]->[6] : $clonehit->{1}->[$i]->[7];
         my $end1   =$strand1 > 0 ? $clonehit->{1}->[$i]->[7] : $clonehit->{1}->[$i]->[6];
         #print "Read1:$chr1\t$strand1\t$start1\t$end1\n";
         ## bes2
         my $chr2=$clonehit->{0}->[$j]->[4];
         my $strand2=$clonehit->{0}->[$j]->[7]-$clonehit->{0}->[$j]->[6] > 0 ? 1 : 0;
         my $start2 =$strand2 > 0 ? $clonehit->{0}->[$j]->[6] : $clonehit->{0}->[$j]->[7];
         my $end2   =$strand2 > 0 ? $clonehit->{0}->[$j]->[7] : $clonehit->{0}->[$j]->[6];
         #print "Read2:$chr2\t$strand2\t$start2\t$end2\n"; 
         ## compare
         my @pos=sort {$a <=> $b} ($start1,$end1,$start2,$end2);
         my $same_chr=$chr1 eq $chr2;
         my $opp_strand=$strand1 ne $strand2;
         my $distance=$same_chr ? abs ($pos[-1]-$pos[0]) : 0;
         if ($same_chr and $opp_strand and $distance >= $min and $distance <= $max){
            push (@hit,[1,$clonehit->{1}->[$i],$clonehit->{0}->[$j]]);  
         }elsif($same_chr and $opp_strand and $distance < $min){
            push (@hit,[2,$clonehit->{1}->[$i],$clonehit->{0}->[$j]]);
         }elsif($same_chr and $opp_strand and $distance > $max){
            push (@hit,[3,$clonehit->{1}->[$i],$clonehit->{0}->[$j]]);
         }elsif($same_chr and !$opp_strand and $distance >= $min and $distance <= $max){
            push (@hit,[4,$clonehit->{1}->[$i],$clonehit->{0}->[$j]]);
         }elsif($same_chr and !$opp_strand and $distance < $min){
            push (@hit,[5,$clonehit->{1}->[$i],$clonehit->{0}->[$j]]);
         }elsif($same_chr and !$opp_strand and $distance > $max){
            push (@hit,[6,$clonehit->{1}->[$i],$clonehit->{0}->[$j]]);
         }elsif($same_chr and $opp_strand and $distance > $maxlimit){
            push (@hit,[7,$clonehit->{1}->[$i],$clonehit->{0}->[$j]]);
         }elsif($same_chr and !$opp_strand and $distance > $maxlimit){
            push (@hit,[8,$clonehit->{1}->[$i],$clonehit->{0}->[$j]]);
         }elsif(!$same_chr){
            push (@hit,[9,$clonehit->{1}->[$i],$clonehit->{0}->[$j]]);
         }
      }
   }
   my @tophit; ### record the highest rank pairs for this clone
   if (@hit == 1){
      push (@tophit,$hit[0]);
   }else{
      my $highest_rank;
      foreach my $h (sort {$a->[0] <=> $b->[0]} @hit){
         if ($highest_rank){
            push (@tophit,$h) if ($h->[0] eq $highest_rank);
         }else{
            push (@tophit,$h);
            $highest_rank=$h->[0]
         }
      }
   }
   ### select the highest score pair in highest rank for final bes pair of this clone.
   my $max_score;
   foreach my $h (@tophit){
      my $tscore =$h->[1]->[12]+$h->[2]->[12]; ## add score of two reads of clone
      $max_score=$tscore > $max_score ? $tscore : $max_score;
   } 
   ### select the lengest pair in highest score pairs
   my $max_length;
   foreach my $h (@tophit){
      my $tlength=$h->[1]->[11]+$h->[2]->[11]; ## add align length of two reads of clone
      my $tscore =$h->[1]->[12]+$h->[2]->[12]; ## add score of two reads of clone
      if ($tscore == $max_score){
         $max_length=$tlength > $max_length ? $tlength : $max_length;
      }
   }
   ### write done the highest score and length pair for this clone
   foreach my $h (@tophit){
      my $tlength=$h->[1]->[11]+$h->[2]->[11];
      my $tscore =$h->[1]->[12]+$h->[2]->[12];
      if ($tscore == $max_score and $tlength == $max_length){
         my $fread=join("\t",@{$h->[1]});
         my $rread=join("\t",@{$h->[2]});
         my $crank=$h->[0];
         #print "$crank\n$fread\n$rread\n";
         $return=[$crank,"$fread\n$rread"];
      }
   }
}
return $return;
}




sub readfpcinf
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $name=$1 if ($unit[0]=~/O\. (\w+)/);
    $hash{$name}=[@unit];
}
close IN;
return \%hash;
}

####
sub statblast
{
my ($file)=@_;
my %clone;
open IN, "$file" or die "$!";
while(my $line=<IN>){
    chomp $line;
    next if ($line=~/^$/);
    my @unit=split("\t",$line);
    my $name=$unit[0];
    if ($name=~/\w+\_+\w+?(a\d+\D+\d+)\.(\w+)/){
       my $clone=$1;
       my $direction=$2;
       #print "$clone\t$direction\n";
       if ($direction eq "r"){
           push @{$clone{$clone}{0}},[@unit];
       }elsif($direction eq "f"){
           push @{$clone{$clone}{1}},[@unit];
       }
    }
}
close IN;
return \%clone;
}



####

