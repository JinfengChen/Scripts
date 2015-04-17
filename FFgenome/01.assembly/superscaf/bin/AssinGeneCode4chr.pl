#!/usr/bin/perl
use Getopt::Long;
use warnings;
#use strict;
GetOptions (\%opt,"gff:s","project:s","help");


my $help=<<USAGE;
Assign gene code to gene and transcripts according to chromosome postion.
We start from ObXXg10010 for each chromosome and increace by 10 every gene.
The transcript id is ObXXg10010.1 for each gene since we have only one transcript for each gene.
perl $0 --gff
--gff
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my ($gene,$genetype,$evidence,$pos,$gene2trans)=grameneGFF($opt{gff});
my $filtergff=$1 if ($opt{gff}=~/(.*)\.gff/);
##sort by position on chromosome
##
my %count;
my $chrcode;
open OUT, ">$filtergff.v1.gff" or die "$!";
foreach my $chr (sort {$a cmp $b} keys %$pos){
  #print "$chr\n";
  my $chrcode=$chr;
  $chrcode=~s/\D+//g;
  $chrcode="Ob$chrcode";
  my @chrgene=@{$pos->{$chr}};
  my @sortchrgene=sort {$a->[1] <=> $b->[1]} @chrgene;
  my $genecode="10000";
  foreach my $g (@sortchrgene){
    $g=$g->[0];
    $genecode+=10;
    my $tempgcode=$chrcode."g".$genecode;
    $count{$genetype->{$g}->[0]}++;
    my @tempgene=split("\t",$genetype->{$g}->[1]);
    $tempgene[8]="ID=$tempgcode;Description=$genetype->{$g}->[0];";
    $tempgene[1]=$opt{project} ? $opt{project} : $tempgene[1];
    my $tempgeneline=join("\t",@tempgene);
    print OUT "$tempgeneline\n";
    #print "$g\n";
    my $tcount=0;
    foreach my $t (keys %{$gene->{$g}}){
       my $len;
       my @mrna;
       my $ref;
       my $strand;
       my $project;
       $tcount++;
       $temptcode=$tempgcode.".".$tcount;
       print "$g\t$t\t$tempgcode\t$temptcode\n";
       foreach my $e (@{$gene->{$g}->{$t}}){
          my @unit=split("\t",$e);
          $ref=$unit[0];
          $project=$opt{project} ? $opt{project} : $unit[1];
          $unit[1]=$project;
          $strand=$unit[6];
          $len+=$unit[4]-$unit[3]+1;
          $unit[8]="Parent=$temptcode;";
          push (@mrna,[@unit]);
       }
       my @sortmrna=sort {$a->[3] <=> $b->[3]} @mrna;
       my $start=$sortmrna[0][3];
       my $end  =$sortmrna[$#sortmrna][4];
       #$len=$end-$start+1;
       my $anno ="ID=$temptcode;"."Parent=$tempgcode;"."Evidence=$evidence->{$t}->[1]";
       unshift (@sortmrna,[$ref,$project,"mRNA",$start,$end,".",$strand,".",$anno]);
       #######write transcript to file
       foreach my $templine (@sortmrna){
          my $line=join("\t",@$templine);
          print OUT "$line\n";
       }
    }
 }
}
close OUT;


####write into file
open OUT, ">$filtergff.summary" or die "$!";
my $totalgene;
foreach(keys %count){
    $totalgene+=$count{$_};
    print OUT "$_\t$count{$_}\n";
}
print OUT "$totalgene\n";
close OUT;
##################

sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    #print "$unit[0]\n";
    $hash{$unit[0]}=1;
}
close IN;
return \%hash;
}

sub grameneGFF
{
my ($gff)=@_;
####read gff formated from gramene gff and store inf as gene->mrna->CDS
my %gene;
my %genetype;
my %evidence;
my %pos;
my %gene2trans;
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/ or $_=~/Component/);
    my $record=$_;
    my @unit=split("\t",$_);
    if ($unit[2]=~/gene/){
        if ($unit[8]=~/ID=(.*?);Description=(.*)$/){
           my $geneline="$unit[0]\t$unit[1]\tgene\t$unit[3]\t$unit[4]\t$unit[5]\t$unit[6]\t$unit[7]\t$unit[8]";
           $genetype{$1}=[$2,$geneline];
           push (@{$pos{$unit[0]}},[$1,$unit[3],$unit[4],$unit[6]]);
        }
    }elsif($unit[2]=~/mRNA/){
        if ($unit[8]=~/ID=(.*?);Parent=(.*?);Evidence=(.*)$/){
            $evidence{$1}=[$2,$3];
            $gene2trans{$2}=$1;
        }
    }elsif($unit[2] eq "CDS"){
        if ($unit[8]=~/ID=(.*?);Parent=(.*?);*$/){
            push (@{$gene{$evidence{$2}->[0]}->{$2}},$record) if (exists $evidence{$2});
        }elsif($unit[8]=~/Parent=(.*?);*$/){
            push (@{$gene{$evidence{$1}->[0]}->{$1}},$record) if (exists $evidence{$1});
            #print "check gff $1\t$evidence{$1}\n";
        }
    }

}
close IN;
return (\%gene,\%genetype,\%evidence,\%pos,\%gene2trans);
}


