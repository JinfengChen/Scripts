#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"CG:s","CHG:s","CHH:s","project:s","help");


my $help=<<USAGE;
Draw CG,CHG,CHH methylation level on gene/repeat and up/downstream.
perl $0 -CG FF.gene.CG.level.sum -CHG FF.gene.CHG.level.sum -CHH FF.gene.CHH.level.sum -project FF.Gene
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}
my %species;
$species{FF}="Oryza brachyantha";
$species{rice}="Oryza sativa";
my ($name,$type,$title);
if ($opt{project}=~/(\w+)\.(\w+)/){
   $name=$1;
   $type=$2;
   $title=$species{$name};
}else{
   $title=$species{$opt{project}};
   $type="Gene";
}
my $max;
my $hei;
my $bin;
my ($lup,$ldown);
if ($type=~/gene/i){
   $max=0.7;
   $hei=0.2;
   $bin=0.1;
}else{
   $max=1;
   $hei=0.4;
   $bin=0.2;
}
$lup=$max-0.02;
$ldown=$max-0.025;

open DRAW, ">$opt{project}.r" or die "$!"; 
print DRAW<<"END.";
pdf("$opt{project}.pdf",width=8,height=5)
scan("$opt{CG}")->x
scan("$opt{CHG}")->y
scan("$opt{CHH}")->z
par(las=1)
num <- c(1:141)
plot(num[1:70],x[1:70],type="l",col="red",font.main=3,lwd=2,axes = FALSE,xlim=c(0,141),ylim=c(0,$max),main="$title",xlab="",ylab="Methylation level")
points(num[72:141],x[71:140],type="l",col="red",lwd=2)
points(num[1:70],y[1:70],type="l",col="blue",lwd=2)
points(num[72:141],y[71:140],type="l",col="blue",lwd=2)
points(num[1:70],z[1:70],type="l",col="lightblue",lwd=2)
points(num[72:141],z[71:140],type="l",col="lightblue",lwd=2)
axis(1,seq(from=0,to=70,by=10),c("-3","-2","-1","0","1","2","3",""))
axis(1,seq(from=71,to=141,by=10),c("","3","2","1","0","-1","-2","-3"))
axis(2,seq(from=0,to=$max,by=$bin),seq(from=0,to=$max,by=$bin))
mtext("4",line=1,side=1,at=70.5)
mtext("$type (kb)",line=2,side=1,at=70)
rect(0,$ldown,3,$lup,border=NA,col="red")
text(6,$lup,"CG",cex=0.7)
rect(12,$ldown,15,$lup,border=NA,col="blue")
text(19,$lup,"CHG",cex=0.7)
rect(25,$ldown,28,$lup,border=NA,col="lightblue")
text(32,$lup,"CHH",cex=0.7)

dev.off()
END.

close DRAW;

`R --vanilla -q <$opt{project}.r`;

