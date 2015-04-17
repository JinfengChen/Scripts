#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"log:s","maxtime:s","project:s","help");


my $help=<<USAGE;
perl $0 --log ./HEG4/make_log/data/run/test/2013-02-13T07\:36\:42/summary.log
--log: summary.log in make_log
--maxtime: max time for y axis, default 200
--project:
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

$opt{project} ||= "HEG4";
$opt{maxtime} ||= 200;

readlog($opt{log});
drawlog();


=pod
[VAPI] Wed Feb 13 07:36:42 PST 2013 : ValidateAllPathsInputs Starting
[VAPI] Wed Feb 13 07:39:30 PST 2013 : ValidateAllPathsInputs Finished
[ln] Wed Feb 13 07:39:30 PST 2013 : ln Starting
[ln] Wed Feb 13 07:39:30 PST 2013 : ln Finished
[ln-2] Wed Feb 13 07:39:30 PST 2013 : ln-2 Starting
[ln-2] Wed Feb 13 07:39:30 PST 2013 : ln-2 Finished
[ln-3] Wed Feb 13 07:39:30 PST 2013 : ln-3 Starting
[ln-3] Wed Feb 13 07:39:30 PST 2013 : ln-3 Finished
[RDR] Wed Feb 13 07:39:30 PST 2013 : RemoveDodgyReads Starting
[RDR] Wed Feb 13 07:50:22 PST 2013 : RemoveDodgyReads Finished
[SPRD] Wed Feb 13 23:59:44 PST 2013 : SamplePairedReadDistributions Starting
[SPRD] Thu Feb 14 00:03:19 PST 2013 : SamplePairedReadDistributions Finished
=cut


sub drawlog
{

my $step=$opt{maxtime}/10;
my $yadj=$opt{maxtime}/50;
my $cmd =<<R;
pdf("$opt{project}.process.pdf",12,7)
par(mai=c(1.7,1,1,0.5))
read.table("$opt{project}.process.time") -> x
barplot(x[,3],ylim=c(0,$opt{maxtime}),axes=FALSE,border=F,col=c("cornflowerblue"),ylab="Time (min)") -> xx
for (i in 1:length(xx)) { # adj can not take vector, so we use loops to add text
  text(xx[i],-$yadj,labels=x[i,2],cex=0.7,srt=65,adj=c(1,1),xpd=TRUE)
}
axis(1,c(0,max (xx)+0.5),line=0.2,labels=c("",""))
axis(2,seq(0,$opt{maxtime},by=$step),cex=1.2)
R
open OUT, ">$opt{project}.process.R" or die "$!";
     print OUT "$cmd\n";
close OUT;
`cat $opt{project}.process.R | R --vanilla --slave`;
}


sub readlog
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
open OUT, ">$opt{project}.process.time" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/ or $_=~/^\[ln/);
    my ($process,$time,$abbr);
    my @unit=split /\s+/,$_;
    $process=$unit[8];
    $abbr=$unit[0];
    $abbr=~s/\[//g;
    $abbr=~s/\]//g;
    #print "$process\t$abbr\n";
    if ($unit[9] eq "Starting"){
       my @unit1=split /\s+/,<IN>;
       last if ($unit1[9] ne "Finished"); 
       if ($unit[2] eq $unit1[2] and $unit[3] eq $unit1[3]){ # difference in one day, time2 - time1
          $time= (time2second($unit1[4]) - time2second($unit[4]))/60 + 1; # add 1 to conpensate the time when int 
          $time= int $time;
          #print "$time\n";
       }elsif($unit[2] eq $unit1[2] and $unit[3] ne $unit1[3]){ # difference bewtten days, days+time1 left in day + time2
          $time=($unit1[3]-$unit[3]-1)*24*60*60 + (24*60*60-time2second($unit[4])) + time2second($unit1[4]);
          $time=int ($time/60);
       }elsif($unit[2] ne $unit1[2]){ # difference bewtten month, days left + time1 left in day + time2
          my $days=monthdays($unit[2]);
          $time=($days-$unit[3])*60*60 + (24*60*60-time2second($unit[4])) + time2second($unit1[4]);
          $time = int ($time/60);
       }
       print OUT "$abbr\t$process\t$time\n";
    }
    
}
close IN;
close OUT;
}
 
sub monthdays
{
my ($month)=@_;
my %data=(
   "Jan" => 31,
   "Feb" => 28,
   "Mar" => 31,
   "Apr" => 30,
   "May" => 31,
   "Jun" => 30,
   "Jul" => 31,
   "Aug" => 31,
   "Sep" => 30,
   "Oct" => 31,
   "Nev" => 30,
   "Dec" => 31
);
return $data{$month};
}


sub time2second
{
my ($time)=@_;
my @unit=split(":",$time);
my $second=$unit[0]*60*60 + $unit[1]*60 + $unit[2];
return $second;
}


