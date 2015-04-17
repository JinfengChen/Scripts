#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"dir:s","help");


my $help=<<USAGE;
Summary gene number and gain/loss information from ml tree and draw figure in R.
This script need to be run in 4tree directory, which contain fa/tree files.
perl $0 --dir ./
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

& DrawR();
print "Type\tOB\tOS\tSB\n";
my @file=glob("$opt{dir}/*.fa");
foreach my $file (@file){
    if ($file=~/$opt{dir}\/(.*)\.fa$/){
       my %gene;
       my $type=$1;
       my $refid=getfastaid($file);
       my ($ob,$os,$sb);
       foreach my $id (keys %$refid){
          if ($id=~/ob/i){
             $ob++;
          }elsif($id=~/os/i){
             $os++;
          }elsif($id=~/sb/i){
             $sb++;
          }
       }
       print "$type\t$ob\t$os\t$sb\n";
       $gene{OBR}=$ob;
       $gene{LOC}=$os;
       $gene{Sb}=$sb;
       my $tree;
       $tree=$file.".muscle.phy_phyml_tree.txt" unless ($type=~/Unknown/);
       $tree=$file.".hmmalign.phy_phyml_tree.txt" if ($type=~/Unknown/);
       & sumtree($tree,\%gene,$type);
    }
}

open OUT, ">>summary.R" or die "$!";
   print OUT <<"END.";
dev.off()
END.
close OUT;


system ("cat summary.R | R --vanilla --slave");


#########################################
sub sumtree
{
my ($tree,$gene,$type)=@_;
`java -jar /home/jfchen/FFproject/tools/Notung-2.6/Notung-2.6.jar -g $tree -s species.txt --speciestag prefix --reconcile --info`;
my $info=$tree.".reconciled.0.info";
#print "$info\n";
my $flag=0;
my $dupli;
my $losses;
my %gain=(
    OBR => 0,
    LOC => 0,
    Sb  => 0,
    Oryza => 0
);

my %loss=(
    OBR => 0,
    LOC => 0,
    Sb  => 0,
    Oryza => 0
);
open IN, "$info" or die "$!";
while(<IN>){
   chomp $_;
   next if ($_=~/^$/);
   my $line=$_;
   if ($line=~/Duplication         L\. Bound            U\. Bound/){
      <IN>;
      $flag=1;
   }elsif($line=~/Duplications: (\d+)/){
      $dupli=$1;
      $flag=2;
   }elsif($line=~/Species             No\. of Losses/){
      <IN>;
      $flag=3;
   }elsif($line=~/Losses: (\d+)/){
      $losses=$1;
   }elsif($flag==1){
      my @unit=split(" ",$line);
      $gain{$unit[1]}++;
   }elsif($flag==3){
      my @unit=split(" ",$line);
      $loss{$unit[0]}+=$unit[1];
   }
}
close IN;

open OUT ,">$type.sum" or die "$!";
   print OUT "OBR\t$gene->{OBR}\t$gain{OBR}\t$loss{OBR}\n";
   print OUT "LOC\t$gene->{LOC}\t$gain{LOC}\t$loss{LOC}\n";
   print OUT "Sb\t$gene->{Sb}\t$gain{Sb}\t$loss{Sb}\n";
   print OUT "Oryza\t0\t$gain{Oryza}\t$loss{Oryza}\n";
close OUT;

my $sum=$type.".sum";
open OUT, ">>summary.R" or die "$!";
   print OUT <<"END.";

sum <- read.table("$sum")
par(adj=0.5,cex=0.7)
plot(x="",xlim=c(0,550),ylim=c(0,200),axes=FALSE,xlab="",ylab="",main="$type")
draw.tree(20,320,140,40,40,sum)

END.
close OUT;
}

sub getfastaid
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
    $hash{$head}=1;
}
$/="\n";
return \%hash;
}

sub DrawR
{
my ($sum,$type)=@_;
open R, ">summary.R" or die "$!";
print R <<"END.";

draw.tree <- function (left,right,height,low,add,sum){
oryza <- left+160
oryzah <- height+add
oryzal <- height-add
lines(x=c(left,right),y=c(low,low))
lines(x=c(left,left),y=c(low,height))
lines(x=c(left,oryza),y=c(height,height))
lines(x=c(oryza,oryza),y=c(oryzal,oryzah))
lines(x=c(oryza,right),y=c(oryzal,oryzal))
lines(x=c(oryza,right),y=c(oryzah,oryzah))
par(adj=0,font=3,cex=0.7)
textx= right+50
text(textx,oryzah,"O. brachyantha")
text(textx,oryzal,"O. sativa")
text(textx,low, "S. bicolor")

###gene number for each species
genex= right+10
text(genex,oryzah,sum[1,2])
text(genex,oryzal,sum[2,2])
text(genex,low, sum[3,2])

###gene number for each node
grassgene=sum[3,2]-sum[3,3]+sum[3,4]
grassx=left+10
grassy=low+50
text(grassx,grassy,grassgene)
oryzagene=sum[2,2]-sum[2,3]+sum[2,4]
oryzax=oryza+10
oryzay=height
text(oryzax,oryzay,oryzagene)

###gene gain and loss from oryza to grass
oryzagl=left+50
oryzagain=height+20
oryzaloss=height-20
text(oryzagl,oryzagain,paste("+",sum[4,3]))
text(oryzagl,oryzaloss,paste("-",sum[4,4]))

### gain and loss for ob
obx= oryza+50
obygain=oryzah+20
obyloss=oryzah-20
text(obx,obygain,paste("+",sum[1,3]))
text(obx,obyloss,paste("-",sum[1,4]))

### gain and loss for os
osx= oryza+50
osygain=oryzal+20
osyloss=oryzal-20
text(osx,osygain,paste("+",sum[2,3]))
text(osx,osyloss,paste("-",sum[2,4]))

### gain and loss for sb
sbx= oryza-20
sbygain=low+20
sbyloss=low-20
text(sbx,sbygain,paste("+",sum[3,3]))
text(sbx,sbyloss,paste("-",sum[3,4]))
}

pdf("summary.pdf")
par(mfcol=c(3,2))
END.
close R;
}


