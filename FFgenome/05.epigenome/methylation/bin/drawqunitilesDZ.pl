#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"max:s","data:s","group:s","type:s","project:s","help");


my $help=<<USAGE;
Draw methylation level for qunitiles grouped by expression level.
--data:  where the bin file is, the results of drawlevel.pl
--group: group file
--max:   max value of y axis
--type:  methylation type
--project: project name
perl $0 -data ./FF -group ../input/Expression/FF.gene.group -type CG -project FF
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

my $bin="$opt{data}/$opt{project}.gene.$opt{type}.level.bin";
my $refbin=parsebin($bin);

my @sumfile;
open IN, "$opt{group}" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    my $rank=shift @unit;
    my $sum =getsum(\@unit,$refbin);
    my $fname="$opt{data}/$opt{project}.$rank.$opt{type}.level.sum";
    push (@sumfile,$fname);
    & write2file($sum,$fname);
}
close IN;

& draw(@sumfile);
########################################
sub draw
{
my (@sumfile)=@_;
my %species=(
     "FF" => "Oryza brachyantha",
     "rice" => "Oryza sativa"
);
my $title=$species{$opt{project}};
my $type="Gene";

my $max=$opt{max};
my $hei;
my $bin;
my ($lup,$ldown);

$max||=0.7;
$hei=0.2;
$bin=0.1;
$lup=$max-0.02;
$ldown=$max-0.025;

open DRAW, ">$opt{data}/$opt{project}.qunitiles.$opt{type}.r" or die "$!";
print DRAW<<"END.";
pdf("$opt{data}/$opt{project}.qunitiles.$opt{type}.pdf",width=8,height=5)
scan("$sumfile[0]")->v
scan("$sumfile[1]")->w
scan("$sumfile[2]")->x
scan("$sumfile[3]")->y
scan("$sumfile[4]")->z
par(las=1)
num <- c(1:141)
plot(num[1:70],v[1:70],type="l",col=1,font.main=3,lwd=2,axes = FALSE,xlim=c(0,141),ylim=c(0,$max),main="$title",xlab="",ylab="methylation level ($opt{type})")
points(num[72:141],v[71:140],type="l",col=1,lwd=2)
points(num[1:70],w[1:70],type="l",col=2,lwd=2)
points(num[72:141],w[71:140],type="l",col=2,lwd=2)
points(num[1:70],x[1:70],type="l",col=3,lwd=2)
points(num[72:141],x[71:140],type="l",col=3,lwd=2);
points(num[1:70],y[1:70],type="l",col=4,lwd=2);
points(num[72:141],y[71:140],type="l",col=4,lwd=2);
points(num[1:70],z[1:70],type="l",col=5,lwd=2);
points(num[72:141],z[71:140],type="l",col=5,lwd=2);

axis(1,seq(from=0,to=70,by=10),c("-3","-2","-1","0","1","2","3",""))
axis(1,seq(from=71,to=141,by=10),c("","3","2","1","0","-1","-2","-3"))
axis(2,seq(from=0,to=$max,by=$bin),seq(from=0,to=$max,by=$bin))
mtext("4",line=1,side=1,at=70.5)
mtext("$type (kb)",line=2,side=1,at=70)

rect(0,$ldown,3,$lup,border=NA,col=1)
text(6,$lup,"1st")
rect(10,$ldown,13,$lup,border=NA,col=2)
text(16,$lup,"2nd")
rect(20,$ldown,23,$lup,border=NA,col=3)
text(26,$lup,"3d")
rect(30,$ldown,33,$lup,border=NA,col=4)
text(36,$lup,"4th")
rect(40,$ldown,43,$lup,border=NA,col=5)
text(46,$lup,"5th")

dev.off()
END.

close DRAW;

`R --vanilla -q < $opt{data}/$opt{project}.qunitiles.$opt{type}.r`;
}

#####
sub getsum
{
my ($genelist,$refbin)=@_;

my @binsum;
my @genesum;
foreach(@$genelist){ 
  my $level=$refbin->{$_};
  for(my $i=0;$i<@$level;$i++){
    unless ($level->[$i] eq "NA"){
            $binsum[$i]+=$level->[$i];
            $genenum[$i]++;
    }
  }
}

my @mean;
for(my $i=0;$i<@binsum;$i++){
   my $mean=$binsum[$i]/$genenum[$i];
   push (@mean,$mean);
}
my $line=join("\t",@mean);
return $line;
}


sub parsebin
{
my ($bin)=@_;
my %hash;
open IN, "$bin" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_ eq "");
    my @unit=split("\t",$_);
    my $gene=shift @unit;
    $hash{$gene}=\@unit;
}
close IN;
return \%hash;
}



sub write2file
{
my ($line,$file)=@_;
open TEMP, ">$file" or die "$!";
    print TEMP "$line\n";
close TEMP;
}

