#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"bed:s","feature:s","methy:s","type:s","project:s","help");


my $help=<<USAGE;
Calculate the methylation level for gene and UTR region by 100 bp interval (5% length for gene body).
--bed:  dir of bed file
--feature: gene/repeat/RT/DNA/COPIA/GYPSY/CACTA/MITE/MUDR
--type:  CG,CHG,CHH
--methy: dir of methylation bed file for CG,CHG,CHH
--project: project name
perl $0 -bed FF_gene_bed -feature gene -type CG -methy FF_methylation_bed --project FF
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

my @genebed=glob("$opt{bed}/*.$opt{feature}.bed");
my @genenum;
my @binsum;
open OUT, ">$opt{project}.$opt{feature}.$opt{type}.level.bin" or die "$!";
open SUM, ">$opt{project}.$opt{feature}.$opt{type}.level.sum" or die "$!";
open PART, ">$opt{project}.$opt{feature}.$opt{type}.level.part" or die "$!";
print PART "ID\tUpC\tUpmC\tBodyC\tBodymC\tDownC\tDownmC\n";
foreach(@genebed){
   print "$_\n";
   my $tempbed=$_;
   my $chr;
   if ($tempbed=~/(\w+)\.$opt{feature}\.bed/){
      $chr=$1;
   }
   my $methbed="$opt{methy}/$chr.$opt{type}.bed";
   next unless (-f $methbed);
   `/home/jfchen/FFproject/tools/BEDTools/bin/windowBed -a $methbed -b $tempbed -w 3000 > chr.win`;
   my ($tempwin,$tempinf)=splitwin("chr.win");
   foreach(keys %$tempwin){
      print OUT "$_\t";
      my $gene=$_;
      my $start=$tempinf->{$_}->[0];
      my $end  =$tempinf->{$_}->[1];
      my $strand=$tempinf->{$_}->[2];
      & tempwin($tempwin->{$_});
      #print "$start\t$end\t$strand\n";
      my ($level,$part)=binlevel("TMP.win",$gene,$start,$end,$strand);
      for(my $i=0;$i<@$level;$i++){
         unless ($level->[$i] eq "NA"){
            $binsum[$i]+=$level->[$i];
            $genenum[$i]++;
         }
         print OUT "$level->[$i]\t";
      }
      print OUT "\n";
      print PART "$gene\t$part->[0]\t$part->[1]\t$part->[2]\t$part->[3]\t$part->[4]\t$part->[5]\n";
   }
}
close OUT;
my @mean;
for(my $i=0;$i<@binsum;$i++){
   my $mean;
   if ($genenum[$i] > 0){
      $mean=$binsum[$i]/$genenum[$i];
   }else{
      $mean="NA";
   }
   push (@mean,$mean);
}
my $line=join("\t",@mean);
print SUM "$line\n";
close SUM;
`rm TMP.win chr.win`;


#######################
sub tempwin
{
my ($line)=@_;
open TEMP, ">TMP.win" or die "$!";
    print TEMP "$line\n";
close TEMP;
}
######
sub splitwin
{
my ($chrwin)=@_;
my %win;
my %inf;
open WIN, "$chrwin" or die "$!";
while(<WIN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    if (exists $win{$unit[10]}){
        $win{$unit[10]}.="$unit[0]\t$unit[1]\t$unit[2]\t$unit[3]\t$unit[4]\t$unit[5]\n";
    }else{
        $inf{$unit[10]}=[$unit[7],$unit[8],$unit[11]];
        $win{$unit[10]}="$unit[0]\t$unit[1]\t$unit[2]\t$unit[3]\t$unit[4]\t$unit[5]\n";
    }
}
close WIN;
return (\%win,\%inf);
}


######
sub binlevel
{
my ($win,$gene,$start,$end,$strand)=@_;
my (@up,@upC,@upmC);
my (@down,@downC,@downmC);
my (@body,@body5C,@body5mC,@body3C,@body3mC);
my @part=(0,0,0,0,0,0);
open WIN, "$win" or die "$!"; 
while(<WIN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    if ($unit[1] < $start){
        my $index=int (abs(($unit[1]-$start)/100));
        $upC[29-$index]++;
        $part[0]++;
        if ($unit[3] > $unit[4]){
           $upmC[29-$index]++;
           $part[1]++;
        }
    }elsif($unit[1] > $end){
        my $index=int (($unit[1]-$end)/100);
        $downC[$index]++;
        $part[4]++;
        if ($unit[3] > $unit[4]){
            $downmC[$index]++;
            $part[5]++;
        }
    }else{
        my $index1=int (($unit[1]-$start)/100);
        $body5C[$index1]++;
        $part[2]++;
        if ($unit[3] > $unit[4]){
           $body5mC[$index1]++;
           $part[3]++;
        }
        my $index2=int (($unit[1]-$end)/100);
        $body3C[39-$index2]++;
        if ($unit[3] > $unit[4]){
           $body3mC[39-$index2]++;
        }
    }
}
close WIN;
my @level;
##up
for(my $i=0;$i<30;$i++){
    my $level;
    if ($upmC[$i] > 0){
       $level=$upmC[$i]/$upC[$i];
    }elsif($upC[$i] > 0){ ### there are C in this region
       $level=0;
    }else{ ### there is no C in this region, so we do not the methylation level
       $level="NA";
    }
    #print "$i\t$upmC[$i]\t$upC[$i]\t$level\n";
    push (@level,$level);
}
### align to five primer
for(my $i=0;$i<40;$i++){
    my $level;
    if ($body5mC[$i] > 0){
       $level=$body5mC[$i]/$body5C[$i];
    }elsif($body5C[$i] > 0){
       $level=0;
    }else{
       $level="NA";
    }
    #print "$i\t$bodymC[$i]\t$bodyC[$i]\t$level\n";
    push (@level,$level);
}
### align to three primer
for(my $i=0;$i<40;$i++){
    my $level;
    if ($body3mC[$i] > 0){
       $level=$body3mC[$i]/$body3C[$i];
    }elsif($body3C[$i] > 0){
       $level=0;
    }else{
       $level="NA";
    }
    #print "$i\t$bodymC[$i]\t$bodyC[$i]\t$level\n";
    push (@level,$level);
}
###down
for(my $i=0;$i<30;$i++){
    my $level;
    if ($downmC[$i] > 0){
       $level=$downmC[$i]/$downC[$i];
    }elsif($downC[$i] > 0){
       $level=0;
    }else{
       $level="NA";
    }
    #print "$i\t$downmC[$i]\t$downC[$i]\t$level\n";
    push (@level,$level);
}

if ($strand eq "-"){
   my @temp=reverse @level;
   undef @level;
   push (@level,@temp);
   my @temp2;
   push (@temp2,@part);
   undef @part;
   $part[0]=$temp2[4];
   $part[1]=$temp2[5];
   $part[2]=$temp2[2];
   $part[3]=$temp2[3];
   $part[4]=$temp2[0];
   $part[5]=$temp2[1];
}
return (\@level,\@part);
}



