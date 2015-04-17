#!/usr/bin/perl
use Getopt::Long;
use strict;
#use warnings;

our %opt;
GetOptions (\%opt,"fasta:s","refer:s","type:s","project:s","help");


my $help=<<USAGE;
Find LRR domain in Pkinase, so the output should be RLK-LRR
perl $0 --fasta Pkinase.OS.fa --project rice --refer rice.RLK-LRR.fa --type rice.RLK-LRR.type
--refer: if set refer, the script will use blast to get sequence
--type: type of refer gene
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

unless ($opt{refer}){
   my $refkinase=findkinase();
   my $reflrr=findLRR();
   my $reflrr2=findLRR2();
   & classify($opt{fasta},$refkinase,$reflrr,$reflrr2);
}else{
   & findRLKLRRbyBLAST();
}
######################################

sub classify
{
my ($fasta,$kinase,$lrr,$lrr2)=@_;
$/="\>";
my $sum;
my %hash;
open OUT, ">$opt{project}.RLK-LRR.fa" or die "$!";
open IN,"$fasta" or die "$!";
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
    next unless (exists $kinase->{$head});
    if (exists $lrr->{$head} or exists $lrr2->{$head}){
       $sum++;
       print OUT ">$head\n$seq\n";
    }
}
print "Total RLK-LRR gene: $sum\n";
$/="\n";
close IN;
close OUT;
`perl /home/jfchen/FFproject/tools/bin/fastaDeal.pl --attr id $opt{project}.RLK-LRR.fa > $opt{project}.RLK-LRR.ID`;
}

##need hmmer2table.pl in this directory
sub findLRR
{
my %hash;
`hmmpfam ../../input/Pfam/LRR.smart.hmm $opt{fasta} > $opt{project}.LRR.hmmpfam`;
`perl hmmer2table.pl --hmmer $opt{project}.LRR.hmmpfam > $opt{project}.LRR.hmmpfam.table`;
open IN, "$opt{project}.LRR.hmmpfam.table" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    next if ($unit[7] > 0.01);
    #print "$unit[0]\t$unit[7]\n";
    $hash{$unit[0]}= $unit[7] < $hash{$unit[0]} ? $unit[7] : $hash{$unit[0]};
}
close IN;
return \%hash;
}

sub findLRR2
{
my %hash;
`hmmsearch --tblout $opt{project}.LRR.hmmsearch.tableout ../../input/Pfam/LRR.hmm $opt{fasta} > $opt{project}.LRR.hmmsearch`;
open IN, "$opt{project}.LRR.hmmsearch.tableout" or die "$!";
     while(<IN>){
            chomp $_;
            next if ($_=~/^#/ or $_=~/^$/);
            my @unit=split(" ",$_);
            next if ($unit[7] > 0.01);
            $hash{$unit[0]}=1;
     }
close IN;
return \%hash;
}

sub findkinase
{
my %hash;
`hmmsearch --tblout $opt{project}.Pkinase.hmmsearch.tableout ../../input/Pfam/Pkinase.hmm $opt{fasta} > $opt{project}.Pkinase.hmmsearch`;
open IN, "$opt{project}.Pkinase.hmmsearch.tableout" or die "$!";
     while(<IN>){
            chomp $_;
            next if ($_=~/^#/ or $_=~/^$/);
            my @unit=split(" ",$_);
            next if ($unit[7] > 0.1);
            $hash{$unit[0]}=1;
     }
close IN;
return \%hash;
}

sub findRLKLRRbyBLAST
{
my %final;
`formatdb -i $opt{fasta}`;
`blastall -p blastp -d $opt{fasta} -i $opt{refer} -o $opt{project}.RLKLRR.blast -e 1e-5 > blast.log 2> blast.log2`;
`/home/jfchen/FFproject/tools/bin/blast_parser.pl $opt{project}.RLKLRR.blast > $opt{project}.RLKLRR.blasttable`;
my $refblast=blastpair("$opt{project}.RLKLRR.blasttable");
my $reftype =readtable($opt{type});
my $refseq  =getfastaseq($opt{fasta});

foreach my $pair (keys %$refblast){
   my @gene=split("\t",$pair);
   #print "$pair\t$refblast->{$pair}\n";
   next if ($refblast->{$pair} < 0.5);
   $final{$reftype->{$gene[0]}}->{$gene[1]}=1;
}

my $total;
foreach my $type (sort keys %final){ 
    open OUT, ">$opt{project}.Blast.$type.fa" or die "$!";
         my $n=keys %{$final{$type}};
         print "$type\t$n\n";
         foreach my $g (keys %{$final{$type}}){
             $total++;
             print OUT ">$g\n$refseq->{$g}\n";
         }
    close OUT;
}
print "Total\t$total\n";
return \%final;
}

sub readtable
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=$unit[1];
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



sub blastpair
{
####read a blasttable file, return a ref of hash contain pair of homologous gene that have e-value <=1e-15,
####and coverage of > 70% for both protein sizes;
my ($blast)=@_;
my %record;
open IN, "$blast" or die "$!";
<IN>;
while(<IN>){
   my @unit=split("\t",$_);
   if ($unit[13] > 1e-15){
       next;
   }
   my $pair=$unit[0]."-".$unit[4];
   if (exists $record{$pair}){
        my $refarray=$record{$pair};
        push (@$refarray,[@unit]);
        $record{$pair}=$refarray;
   }else{
        my @array;
        push (@array,[@unit]);
        $record{$pair}=\@array;
   }
}
close IN;

my %homologous;
foreach (keys %record){
      #print "$_\n";
      #print Dumper($record{$_}),"\n";
      my $refarray=$record{$_};
      my $num=@$refarray;
      #print "$num\n";
      my $hsp=1;
      my $queryid =$$refarray[0][0];
      my $hitid   =$$refarray[0][4];
      my $bitscore=$$refarray[0][12];
      my $identity=$$refarray[0][8];
      my $querylen=$$refarray[0][1];
      my $hitlen  =$$refarray[0][5];
      my $longestq=$$refarray[0][3]-$$refarray[0][2]+1;
      my $longesth=$$refarray[0][7]-$$refarray[0][6]+1;
      my $matchq  =$$refarray[0][3]-$$refarray[0][2]+1;
      my $matchh  =$$refarray[0][7]-$$refarray[0][6]+1;
      if ($num > 1){
          my $hspstartq=$$refarray[0][2];
          my $hspstarth=$$refarray[0][6];
          my $hspendq  =$$refarray[0][3];
          my $hspendh  =$$refarray[0][7];
          for(my $i=1;$i<=$num-1;$i++){
                if ($$refarray[$i][2] < $hspendq or $$refarray[$i][6] < $hspendh){
                      next;
                }else{
                      $hsp++;
                      $hspendq=$$refarray[$i][3];
                      $hspendh=$$refarray[$i][7];
                      $bitscore+=$$refarray[$i][12];
                      $identity+=$$refarray[$i][8];
                      $matchq  +=$$refarray[$i][3]-$$refarray[$i][2]+1;
                      $matchh  +=$$refarray[$i][7]-$$refarray[$i][6]+1;
                }
          }
          $longestq=$hspendq-$hspstartq+1;
          $longesth=$hspendh-$hspstarth+1;
          $identity=$identity/$hsp;
      }
      my $qhspcoverage=$matchq/$querylen;
      my $hhspcoverage=$matchh/$hitlen;
      my $qlencoverage=$longestq/$querylen;
      my $hlencoverage=$longesth/$hitlen;
      #if ($qhspcoverage >= 0.7 and $hhspcoverage >= 0.7){
      #if ($qlencoverage >= 0.7 and $hlencoverage >= 0.7){
      if ($qhspcoverage >= 0.3 and $hhspcoverage >= 0.3 and $qlencoverage >= 0.5 and $hlencoverage >= 0.5){
          my $temp="$queryid\t$hitid";
          #print "$temp\n";
          #next if ($identity < 60);
          $homologous{$temp}=$identity;
      }
}
return \%homologous;
}

