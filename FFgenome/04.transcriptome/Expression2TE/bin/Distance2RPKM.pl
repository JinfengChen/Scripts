#!/usr/bin/perl
use Getopt::Long;
use Data::Dumper;

GetOptions (\%opt,"bed:s","gff:s","rpkm:s","me:s","type:s","project:s","help");


my $help=<<USAGE;
Relation of distance from nearest TE and Expression level for gene in one species.
-bed: BED format gene annotation file, generate by GFF2BED.pl from GFF3 format gene annotation.
-gff: BED format TE annotation file, should be contain only TE class that interested to calculate correlation.
-rpkm: gene expression data.
Gene	RPKM
OBR_GLEAN_10017382	64.3014830721204
OBR_GLEAN_10011719	78.1501432555481
OBR_GLEAN_10007059	24.7132214455825
-project: project name.
-me: methylation status(id,status)
-type: all/upstream/downstream closest feature
Run: perl $0 -bed FF.mRNA.bed -gff ../input/FF.repeat.gff -rpkm ../input/FF.shoot.rpkm -me ../input/FF.repeat.gff.Me.status -type all -project FF
 
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my %hash; ## change this hash to meth and unmeth to store distance and expression
my %me;
my %unme;
my @control;
my $control;
my $max=1000;
my $bin=$max/10;
my $controlindex=500/$bin;
my $cutoff=0.1;
my $BEDtools="/home/jfchen/FFproject/tools/BEDTools/bin";
my $refrpkm=rpkmexpr($opt{rpkm});
my $refme  =mestatus($opt{me});

my $refBED;
if ($opt{type}=~/all/i){
    system("$BEDtools/closestBed -a $opt{bed} -b $opt{gff} > $opt{project}.closestBED");
    $refBED=closestBED("$opt{project}.closestBED");
    `rm *.closestBED`;
}elsif($opt{type}=~/up|5/){
    print "Get upstream cloest feature\n";
    $refBED=upstream($opt{bed},$opt{gff});
}elsif($opt{type}=~/down|3/){
    print "Get downstream cloest feature\n";
    $refBED=downstream($opt{bed},$opt{gff});
}
open OUT, ">$opt{project}.gene.closestTE.inf" or die "$!";
print OUT "Gene\tRPKM\tRepeat\tDistance\tMethylation\n";
foreach(keys %$refBED){
    my $rpkm0= $refrpkm->{$_}->[0] ? $refrpkm->{$_}->[0] : "NA";
    my $rpkm= $refrpkm->{$_}->[1] ? $refrpkm->{$_}->[1] : "NA";
    print OUT "$_\t$rpkm0\t$refBED->{$_}->[0]\t$refBED->{$_}->[1]\t$refme->{$refBED->{$_}->[0]}\n";
    next if ($rpkm eq "NA");
    my $index;
    if ($refBED->{$_}->[1] > 0){
       $index=int ($refBED->{$_}->[1]/$bin) + 1;
    }else{
       $index=0;
    }
    my $dist=$index*$bin;
    next if ($dist > $max);
    push (@{$hash{$index}},$rpkm); ## test if methylation or not and push to meth or unmeth
    if ($refme->{$refBED->{$_}->[0]} >= $cutoff){
         push (@{$me{$index}},$rpkm);
    }elsif($refme->{$refBED->{$_}->[0]} < $cutoff and $refme->{$refBED->{$_}->[0]} ne "NA"){
         push (@{$unme{$index}},$rpkm);
    }
    if ($index > $controlindex){
       push (@control,$rpkm);
    }
}
close OUT;
$control=mean(\@control);

open OUT, ">$opt{project}.all.4r" or die "$!";
my @pos;
my @all;
foreach(sort {$a <=> $b} keys %hash){
    my $array=$hash{$_};
    my $temp=@$array;
    my $dis =$_*$bin;
    push (@pos,$dis);
    push (@all,$temp);
    my $dist =$_*$bin;
    foreach(@$array){
        print OUT "$dist\t$_\n";
    }
}
close OUT;

my @me;
open OUT, ">$opt{project}.me.4r" or die "$!";
foreach(sort {$a <=> $b} keys %me){
    my $array=$me{$_};
    my $dist =$_*$bin;
    my $temp=@$array;
    push (@me,$temp);
    foreach(@$array){
        print OUT "$dist\t$_\n";
    }
}
close OUT;

my @unme;
open OUT, ">$opt{project}.unme.4r" or die "$!";
foreach(sort {$a <=> $b} keys %unme){
    my $array=$unme{$_};
    my $temp=@$array;
    push (@unme,$temp);
    my $dist =$_*$bin;
    foreach(@$array){
        print OUT "$dist\t$_\n";
    }
}
close OUT;


open SUM, ">$opt{project}.summary" or die "$!";
print SUM "Distance\tAll\tMe\tUnMe\n";
for(my $i=0;$i<@all;$i++){
   print SUM "$pos[$i]\t$all[$i]\t$me[$i]\t$unme[$i]\n";
}
close SUM;


open OUT, ">$opt{project}.r" or die "$!"; 
print OUT <<"END.";
read.table("$opt{project}.all.4r") -> x
read.table("$opt{project}.me.4r") -> y
read.table("$opt{project}.unme.4r") -> z
pdf("$opt{project}.pdf")
library("gplots")
plotmeans(y[,2]~y[,1],ylim=c(-1,1),n.label=FALSE,barcol="black",ylab="Expression (RPKM)",xlab="Distance to gene (bp)")
par(new=T)
plotmeans(z[,2]~z[,1],ylim=c(-1,1),n.label=FALSE,col="red",barcol="red",xlab="",ylab="",new=TRUE);
abline(h=$control)
dev.off()
END.
close OUT;

system ("cat $opt{project}.r | R --vanilla --slave");
################################################
sub closestBED
{
#### get the distance from gene to nearest TE 
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    my $distance;
    my $id;
    if ($unit[9] > $unit[2]){
       $distance=$unit[9]-$unit[2]+1;
    }elsif($unit[10] < $unit[1]){
       $distance=$unit[1]-$unit[10]+1;
    }else{
       $distance=0;
    }
    if ($unit[14]=~/ID=(.*?);/){
       $id=$1;
    }
    $hash{$unit[3]}=[$id,$distance];
}
close IN;
return \%hash;
}

sub mestatus
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    $hash{$unit[0]}=$unit[1];
}
close IN;
return \%hash;
}

sub rpkmexpr
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    $hash{$unit[0]}=[$unit[1],$unit[2]];
}
close IN;
return \%hash;
}

=pod
sub rpkmexpr
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split(" ",$_);
    if ($unit[1]=~/\-/ or $unit[1] == 0){
       $hash{$unit[0]}=[-2,0.01];
       #$hash{$unit[0]}=["NA","NA"];
    }else{
       my $temp=log10($unit[1]);
       $hash{$unit[0]}=[$temp,$unit[1]];
    }
}
close IN;
return \%hash;
}
=cut

sub mean
{
my ($num)=@_;
my $loop=0;
my $total;
foreach  (@$num) {
        next if ($_ eq "NA");
        $total+=$_;
        $loop++;
}
$number=@$num;
$mean=$total/$loop;
$median=$num->[int $number/2];
return $mean;
}

sub log10 {
    my ($n) = shift;
    return log($n)/log(10);
}
##closest upstream feature
sub upstream
{
my ($bed,$gff)=@_;
my %hash;
my $refbed=bed($bed);
my ($refstart,$refend)=gff($gff);
print "Gene\tStrand\tID\tStart\tCloestend\tIndex\tDistance\n";
foreach(keys %$refbed){
   my $gene=$_;
   my $strand=$refbed->{$gene}->[0][3];
   print "$gene\t$strand\t"; 
   if ($refbed->{$gene}->[0][3] =~/\+/ ){
      my $start=$refbed->{$gene}->[0][1];
      my $teend=$refend->{$refbed->{$gene}->[0][0]};
      my @end=keys %$teend;
      push (@end,$start);
      my @array = sort { $a <=> $b} @end;
      my $index = bsearch(\@array,$start);
      next if ($index == 0 or $index == $#array); ## no match, may at the chromosome end
      #print "$array[0]\t$array[$#array]\n";
      #my ($index) = grep { $array[$_] == $start } 0..$#array;
      my $closestend = $array[$index-1];
      my $teid =$teend->{$closestend};
      my $distance =$start-$closestend+1;
      print "$teid\t$start\t$closestend\t$index\t$distance\n";
      $hash{$gene}=[$teid,$distance]; 
   }else{
      my $start=$refbed->{$gene}->[0][2];
      my $testart=$refstart->{$refbed->{$gene}->[0][0]};
      my @start=keys %$testart;
      push (@start,$start);
      my @array = sort { $a <=> $b} @start;
      my $index = bsearch(\@array,$start);
      next if ($index == 0 or $index == $#array); ## no match, may at the chromosome end
      #my ($index) = grep { $array[$_] == $start } 0..$#array;
      my $closeststart = $array[$index+1];
      my $teid =$testart->{$closeststart};
      my $distance =$closeststart-$start+1;
      print "$teid\t$start\t$closeststart\t$index\t$distance\n";
      $hash{$gene}=[$teid,$distance];
   }
}
return \%hash;
}
##closest downstream feature
sub downstream
{
my ($bed,$gff)=@_;
my %hash;
my $refbed=bed($bed);
my ($refstart,$refend)=gff($gff);
print "Gene\tStrand\tID\tStart\tCloestend\tIndex\tDistance\n";
foreach(keys %$refbed){
   my $gene=$_;
   my $strand=$refbed->{$gene}->[0][3];
   print "$gene\t$strand\t";
   unless ($refbed->{$gene}->[0][3] =~/\+/ ){
      my $start=$refbed->{$gene}->[0][1];
      my $teend=$refend->{$refbed->{$gene}->[0][0]};
      my @end=keys %$teend;
      push (@end,$start);
      my @array = sort { $a <=> $b} @end;
      my $index = bsearch(\@array,$start);
      next if ($index == 0 or $index == $#array); ## no match, may at the chromosome end
      #my ($index) = grep { $array[$_] == $start } 0..$#array;
      my $closestend = $array[$index-1];
      my $teid =$teend->{$closestend};
      my $distance =$start-$closestend+1;
      print "$teid\t$start\t$closestend\t$index\t$distance\n";
      $hash{$gene}=[$teid,$distance];
   }else{
      my $start=$refbed->{$gene}->[0][2];  ## $start here is end of gene on + strand
      my $testart=$refstart->{$refbed->{$gene}->[0][0]};
      my @start=keys %$testart;
      push (@start,$start);
      my @array = sort { $a <=> $b} @start;
      my $index = bsearch(\@array,$start);
      next if ($index == 0 or $index == $#array); ## no match, may at the chromosome end
      #my ($index) = grep { $array[$_] == $start } 0..$#array;
      my $closeststart = $array[$index+1];
      my $teid =$testart->{$closeststart};
      my $distance =$closeststart-$start+1;
      print "$teid\t$start\t$closeststart\t$index\t$distance\n";
      $hash{$gene}=[$teid,$distance];
   }
}
return \%hash;
}

sub gff
{
my ($gff)=@_;
my %hash;  ##te start
my %hash2; ##te end
my $seq;
my $id;    ##ID for element
open IN, "$gff" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    $seq=$unit[0];
    if ($unit[8]=~/ID=(.*?);/){
       $id=$1;
    }
    unless (exists $hash{$seq}){
       my %testart;
       $testart{$unit[3]}=$id;
       $hash{$seq}=\%testart;
       my %teend;
       $teend{$unit[4]}=$id;
       $hash2{$seq}=\%teend;
    }else{
       my $testart=$hash{$seq};
       $testart->{$unit[3]}=$id;
       $hash{$seq}=$testart;
       my $teend=$hash2{$seq};
       $teend->{$unit[4]}=$id;
       $hash2{$seq}=$teend;
    }

}
close IN;
return (\%hash,\%hash2);
}


sub bed
{
my ($bed)=@_;
my %hash;  ##gene inf
open IN, "$bed" or die "$!";
while (<IN>){
    chomp $_;
    next if ($_=~/^#/);
    my @unit=split("\t",$_);
    push (@{$hash{$unit[3]}},[$unit[0],$unit[1],$unit[2],$unit[5]]);
}
close IN;
return \%hash;
} 

sub bsearch 
{
    my ($array, $word) = @_;
    my $low = 0;
    my $high = @$array - 1;

    while ( $low <= $high ) {
        my $try = int( ($low+$high) / 2 );
        $low  = $try+1, next if $array->[$try] < $word;
        $high = $try-1, next if $array->[$try] > $word;
        return $try;
    }
    return;
}

