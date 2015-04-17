#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"blastm8:s","repeat:s","help");


my $help=<<USAGE;
perl $0 --blastm8 --repeat
--blastm8: name.blast.m8
--repeat: gff of reference repeat annotation
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

`perl /home/jfchen/FFproject/tools/solar/solar/solar.pl -a est2genome1 -f m8 -z $opt{blastm8} > exon.solar`;
`perl /home/jfchen/FFproject/tools/bin/bestAlign.pl exon.solar > exon.best.solar`;
my $position=findbadexon("exon.best.solar");
my $te=parseGFFrepeat($opt{repeat});
& gff2table($opt{repeat},"Trans","repeat.table");
`perl findOverlap.pl gapexon.txt repeat.table > repeat.overlap`;
my ($exon,$repeat)=sumrepeat("repeat.overlap",$te,$position);

########
sub findbadexon
{
my ($file)=@_;
my %position;
open OUT, ">gapexon.txt" or die "$!";
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $exon=$unit[0];
    my $chr =$unit[5];
    $position{$exon}=[$chr,$unit[7],$unit[8]];
    next if $unit[9] < 2;
    my @qryblock=split(";",$unit[11]);
    my @refblock=split(";",$unit[12]);
    my @gap;
    my $flag=0;
    for(my $i=1;$i<@qryblock;$i++){
       my @ref=split(",",$refblock[$i]);
       my @ref0=split(",",$refblock[$i-1]);
       my $reflen=abs($ref[0]-$ref0[1]);
       my @qry=split(",",$qryblock[$i]);
       my @qry0=split(",",$qryblock[$i-1]);
       my $qrylen=$qry[0]-$qry0[1];
       my $start=$ref[0] > $ref0[1] ? $ref0[1] : $ref[0];
       my $end  =$ref[0] > $ref0[1] ? $ref[0] : $ref0[1];
       #print "$exon\t$qrylen\t$reflen\n";
       if ($reflen > 200){
          #push (@gap,"$chr:$start:$end:$reflen:$qrylen");
          push (@gap,"$chr\t$exon\t$start\t$end\t$reflen:$qrylen");
          $flag=1;
       }
    }
    my $hit=join("\n",@gap);
    print OUT "$hit\n" if ($flag == 1);
}
close IN;
close OUT;
return \%position;
}

sub gff2table
{
my ($gff,$feature,$table)=@_;
open IN, "$gff" or die "$!";
open OUT, ">$table.unsort" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    if ($unit[2]=~$feature and ($unit[8] =~/ID=(.*?);/ or $unit[8] =~/Parent=(.*?);/)){
        print OUT "$unit[0]\t$1\t$unit[3]\t$unit[4]\n";
    }
}
close IN;
close OUT;
system ("msort -k 1,n3 $table.unsort > $table");
system ("rm $table.unsort");
}

sub sumrepeat
{
my ($file,$te,$position)=@_;
my %hash;
my %exon;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split(" ",$_);
    next if ($unit[3] == 0);
    my %temp;
    my $TElen;
    for(my $i=4;$i<@unit;$i++){
       my @array=split(",",$unit[$i]);
       $temp{$array[0]}=$array[2];
       $TElen+=$array[2];
    } 
    my @id=sort {$temp{$b} <=> $temp{$a}} keys %temp;
    my $id=$id[0];
    my $len=$temp{$id};
    my $type=$te->{$id};
    print "$unit[0]\t$position->{$unit[0]}->[0]\t$position->{$unit[0]}->[1]\t$position->{$unit[0]}->[2]\t$unit[1]\t$id\t$type\t$len\t$TElen\n" if ($TElen/$unit[1] > 0.7);
    $hash{$type}+=$len;
}
close IN;
return \%hash;
}


sub parseGFFrepeat
{
my ($gff)=@_;
my %hash;  ##hash to store every record by key of Seq_id
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    if ($unit[8]=~/ID=(.*?);.*Class=(.*?);/){
       my $id=$1;
       my $class=$2;
       if ($class=~/LTR/){
          $class="LTR";
       }elsif($class=~/DNA/){
          $class="DNA";
       }else{
          $class="otherTE";
       }
       $hash{$id}=$class
    }
}
close IN;
return \%hash;
}


###########################
sub getblat
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $exon=$unit[9];
    my @size=split(",",$unit[18]);
    my @qrystart=split(",",$unit[19]);
    my @refstart=split(",",$unit[20]);
    next if @size < 2;
    my @gap;
    my $flag=0;
    for(my $i=1;$i<@size;$i++){
       my $qrylen=$qrystart[$i]-($qrystart[$i-1]+$size[$i-1]);
       my $reflen=$refstart[$i]-($refstart[$i-1]+$size[$i-1]);
       if ($reflen > 100 and $qrylen < 20){
          $flag=1;
          my $start=$refstart[$i-1]+$size[$i-1];
          my $end  =$refstart[$i];
          my $chr  =$unit[13];
          push (@gap,"$chr:$start:$end");
       } 
    }
    my $hit=join("\t",@gap);
    print "$exon\t$hit\n" if ($flag == 1);
}
close IN;
return \%hash;
}

