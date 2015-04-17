#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"pfamsum:s","project:s","help");


my $help=<<USAGE;
perl $0 --pfamsum TIGR6.Pfam.summary --project TIGR6
--pfamsum:
--project: project name.
 
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

& barplot($opt{pfamsum},"1","4","0.6","Non Ortholog","$opt{project}.nonorth");
& barplot($opt{pfamsum},"2","4","0.6","Non-synteny Ortholog","$opt{project}.nonsynorth");
& barplot($opt{pfamsum},"3","4","1","Synteny Ortholog","$opt{project}.synorth");

sub barplot
{
my ($file,$col1,$col2,$ylim,$ylab,$title)=@_;
my %hash;
my %marker=(
   "10" => "<10",
   "50" => "<50",
   "100" => "<100",
   "200" => "<200",
   "500" => "<500",
   "1000" => ">=500"
);
open IN, "$file" or die "$!";
<IN>;
while(<IN>){
   chomp $_;
   next if ($_=~/^$/);
   my @unit=split("\t",$_);
   if ($unit[$col2] <= 10){
      my $rate=$unit[$col1]/$unit[$col2];
      push (@{$hash{"10"}},$rate);
   }elsif($unit[$col2] < 50){
      my $rate=$unit[$col1]/$unit[$col2];
      push (@{$hash{"50"}},$rate);
   }elsif($unit[$col2] < 100){
      my $rate=$unit[$col1]/$unit[$col2];
      push (@{$hash{"100"}},$rate);
   }elsif($unit[$col2] < 200){
      my $rate=$unit[$col1]/$unit[$col2];
      push (@{$hash{"200"}},$rate);
   }elsif($unit[$col2] < 500){
      my $rate=$unit[$col1]/$unit[$col2];
      push (@{$hash{"500"}},$rate);
   }else{
      my $rate=$unit[$col1]/$unit[$col2];
      push (@{$hash{"1000"}},$rate);
   }
}
close IN;

open OUT, ">$title.4r" or die "$!";
foreach(sort {$a <=> $b} keys %hash){
    my $size=$marker{$_};
    my ($mean,$se,$number)=ci($hash{$_});
    my $cilow=$mean-$se;
    my $cihigh=$mean+$se;
    print OUT "$mean\t$cilow\t$cihigh\t$size\t$se\n";
}
close OUT;

open OUT, ">$title.r" or die "$!";
print OUT <<"END.";
read.table("$title.4r") -> x
pdf("$title.pdf")

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop("vectors must be same length")
arrows(x,y+upper, x, y-upper, angle=90, code=3, length=length, ...)
}

barx <- barplot(x[,1],ylim=c(0,$ylim),names.arg=x[,4],ylab="$ylab",xlab="Pfam Family Size")
error.bar(barx,x[,1],x[,5]) 
abline(h=0)
dev.off()

END.
close OUT;

system ("cat $title.r | R --vanilla --slave");
}




sub ci
{
my ($num)=@_;
my $loop=0;
my $total=0;
my $add_square;
foreach  (@$num) {
        next if ($_ eq "NA");
        #my $temp=log($_);
        #print "$_\t$temp\n";
        my $temp=$_;
        $total+=$temp;
        $add_square+=$temp*$temp;
        $loop++;
}


my $number=$loop;
return ($num->[0],0,1) if ($number < 2);
my $mean=$total/$number;
my $SD=sqrt( ($add_square-$total*$total/$number)/ ($number-1) );
my $se=1.96*$SD/sqrt($number);
return ($mean,$se,$number);
}


 
