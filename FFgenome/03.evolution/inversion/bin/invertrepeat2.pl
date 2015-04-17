#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"inversion:s","refseq:s","qryseq:s","refgff:s","qrygff:s","flanking:s","help");


my $help=<<USAGE;
perl $0 --inversion --refseq --qryseq 
Analysis flanking sequence of inversion to find invert repeat.
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

$opt{flanking} ||= 10000;
my $refseq=getfastaseq($opt{refseq});
my $qryseq=getfastaseq($opt{qryseq});
my $refgff=parseRepeatpos($opt{refgff});
my $qrygff=parseRepeatpos($opt{qrygff});
=pod
foreach my $c (keys %$refgff){
   print "$c\n";
   my $t=$refgff->{$c};
   foreach my $te (keys %{$refgff->{$c}}){
      print "TE:$te\n";
   }
}
=cut
my $rank;
open OUT, ">IR.list" or die "$!";
open IN, "$opt{inversion}" or die "$!";
while(<IN>){
   chomp $_;
   next if ($_=~/^$/);
   $rank++;
   my @unit=split("\t",$_);
   my ($refs1,$refs2)=invert($refseq->{"chr$unit[0]"},$unit[5],$unit[6],"ref$rank");
   print "$rank\nRef\tchr$unit[0]\t$unit[5]\t$unit[6]\t$refs1\t$refs2\n";
   my $hitref=ssearch($refgff->{"chr$unit[0]"},$refs1,$refs2,"ssearch.out");
   my ($qrys1,$qrys2)=invert($qryseq->{"chr$unit[0]"},$unit[7],$unit[8],"qry$rank");
   print "Qry\tchr$unit[0]\t$unit[7]\t$unit[8]\t$qrys1\t$qrys2\n";
   my $hitqry=ssearch($qrygff->{"chr$unit[0]"},$qrys1,$qrys2,"ssearch.out");
   print OUT "$rank\t$unit[1]\t$hitref\t$hitqry\n";
}
close IN;
close OUT;

#######################
sub invert
{
my ($seq,$temp1,$temp2,$title)=@_;
my $start=$temp1 > $temp2 ? $temp2 : $temp1;
my $end  =$temp1 > $temp2 ? $temp1 : $temp2;
my $s1=$start-$opt{flanking};
my $subseq1=substr($seq,$s1,$opt{flanking});
writefasta("$title.up.fa",$subseq1);
my $s2=$end;
my $subseq2=substr($seq,$s2,$opt{flanking});
writefasta("$title.down.fa",$subseq2);
`ssearch36 $title.up.fa $title.down.fa -m 8 > ssearch.out`;
`rm *.fa`;
#print "S1S2:$s1\t$s2\n";
return ($s1,$s2);
}

sub ssearch
{
my ($gff,$s1,$s2,$file)=@_;
my $hit;
open IN1, "$file" or die "$!";
while(<IN1>){
    chomp $_;
    next if ($_=~/^$/ or $_=~/^#/);
    my @unit=split("\t",$_);
    my $flag1=$unit[7]-$unit[6] > 0 ? "+" : "-";
    my $flag2=$unit[9]-$unit[8] > 0 ? "+" : "-";
    if ($flag1 ne $flag2 or $unit[10] < 0.01){
       my $m1u=$unit[7]-$unit[6] > 0 ? $unit[6]+$s1 : $unit[7]+$s1;
       my $m1d=$m1u+abs($unit[7]-$unit[6]);
       my $overlap1=overlap($gff,$m1u,$m1d,$flag1);
       print "Up Hit:$m1u\t$m1d\t$flag1\n$overlap1\n";
       my $m2u=$unit[9]-$unit[8] > 0 ? $unit[8]+$s2 : $unit[9]+$s2;
       my $m2d=$m2u+abs($unit[9]-$unit[8]);
       my $overlap2=overlap($gff,$m2u,$m2d,$flag2);
       print "Down Hit:$m2u\t$m2d\t$flag2\n$overlap2\n";
       $hit="$overlap1,$overlap2";
       last;
    } 
}
close IN1;
return $hit;
}

sub overlap
{
my ($gff,$mu,$md,$strand)=@_;
my @overlap;
foreach my $te (keys %$gff){
   #print "TE:$te\t$gff->{$te}->[0]\t$gff->{$te}->[1]\n";
   if ($gff->{$te}->[0] > $mu and $gff->{$te}->[0] < $md){
      my $class=$gff->{$te}->[3];
      my $overlapsize= $gff->{$te}->[1] > $md ? $md-$gff->{$te}->[0]+1 : $gff->{$te}->[1]-$gff->{$te}->[0]+1;
      push (@overlap,[$te,$class,$gff->{$te}->[0],$gff->{$te}->[1],$overlapsize]);
   }elsif($gff->{$te}->[1] > $mu and $gff->{$te}->[1] < $md){
      my $class=$gff->{$te}->[3];
      my $overlapsize= $mu > $gff->{$te}->[0] ? $gff->{$te}->[1]-$mu+1 : $gff->{$te}->[1]-$gff->{$te}->[0]+1;
      push (@overlap,[$te,$class,$gff->{$te}->[0],$gff->{$te}->[1],$overlapsize]);
   }
}
my @sort=sort {$b->[4] <=> $a->[4]} @overlap;
my $longest="$sort[0][0]:$sort[0][1]:$sort[0][2]:$sort[0][3]:$sort[0][4]"; ##id:class:start:end:overlap
return $longest;
}


sub writefasta
{
my ($file,$seq)=@_;
my $name=$1 if ($file=~/(.*)\.fa/);
open OUT1, ">$file" or die "$!";
     print OUT1 ">$name\n$seq\n";
close OUT1;
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
close IN;
$/="\n";
return \%hash;
}

##ID=1;Target=FPH184:1..190;Class=DNA/Harbinger;
sub parseRepeatpos
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
my $id;    ##ID for element
my $type;
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[2]=~/Trans/){
        if ($unit[8]=~/ID=(.*?);.*;Class=(.*?);/){
            $id=$1;
            $type=$2;
        }
        $hash{$unit[0]}{$id}=[$unit[3],$unit[4],$unit[6],$type]; 
    }
}
close IN;
return \%hash;
}

