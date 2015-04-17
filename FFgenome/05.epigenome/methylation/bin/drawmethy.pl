#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"CG:s","CHG:s","CHH:s","project:s","help");


my $help=<<USAGE;
Draw CG,CHG,CHH methylation level on gene/repeat and up/downstream.
perl $0 -CG FF.CG.level.sum -CHG FF.CHG.level.sum -CHH FF.CHH.level.sum -project FF.Gene
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
   $max=0.5;
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
pdf("$opt{project}.pdf",width=6,height=4)
scan("$opt{CG}")->x
scan("$opt{CHG}")->y
scan("$opt{CHH}")->z
par(las=1)
plot(x,type="l",col="red",font.main=3,lwd=2,axes = FALSE,xlim=c(0,60),ylim=c(0,$max),main="$title",xlab="",ylab="Methylation level")
points(y,type="l",col="blue",lwd=2)
points(z,type="l",col="lightblue",lwd=2)
axis(1,seq(from=0,to=60,by=10),c(2,1,0,50,100,1,2))
axis(2,seq(from=0,to=$max,by=$bin),seq(from=0,to=$max,by=$bin))
rect(20,-$hei,40,0,col="black")
mtext("Upstream (kb)",line=2,side=1,at=10)
mtext("$type (%)",line=2,side=1,at=30)
mtext("Downstream (kb)",line=2,side=1,at=50)
rect(0,$ldown,3,$lup,border=NA,col="red")
text(5,$lup,"CG")
rect(10,$ldown,13,$lup,border=NA,col="blue")
text(16,$lup,"CHG")
rect(20,$ldown,23,$lup,border=NA,col="lightblue")
text(26,$lup,"CHH")

dev.off()
END.

close DRAW;

`R --vanilla -q <$opt{project}.r`;

