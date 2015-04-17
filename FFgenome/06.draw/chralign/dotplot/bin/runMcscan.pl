#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"mcl","project:s","help");

my $help=<<USAGE;
Run this program in bin dir, cp os_sb.* into ./bin

Run: perl runMcscan.pl -mcl -project ./ob_os or perl runMcscan.pl -mcl -project ob_os

-mcl: if do cluster using mcl, if mcl is not set we only do pairware align using -a.
-project: prefix of blast and gff

USAGE

if ($opt{help} or keys %opt < 1){
   print "$help\n";
   exit();
}


#my $qsub    ="/share/raid12/chenjinfeng/tools/bin/qsub-sge.pl";
#my $mcscan  ="/share/raid12/chenjinfeng/tools/mcscan0.8/mcscan";
#my $mcscan  ="./mcscan"; # this one use gff
my $mcscan  ="./mcscan_last"; # this one do use bed
my $mcl     ="./mcl";
my $mc2     ="./mc2dagalign.pl";
my $dot     ="./bcdotplot.pl";

my $mcscanshell="mcscanshell.sh";
if ($opt{mcl}){
   system "more $opt{project}.blast | $mcl - --abc --abc-neg-log -abc-tf 'mul(0.4343), ceil(200)' -o $opt{project}.mcl > mcl.log 2> mcl.log2 ";
   #open OUT, ">$mcscanshell" or die "$!";
   #    print OUT "$mcscan $opt{project} > log 2> log2 &";
   #close OUT;
   #system "$qsub --resource vf=0.4G $mcscanshell &";
   system "$mcscan -p Os -s 5 $opt{project} > log 2> log2";
}else{
   #open OUT, ">$mcscanshell" or die "$!";
   #    print OUT "$mcscan $opt{project} -a -g -3 > log 2> log2 &";
   #close OUT;
   #system "$qsub --resource vf=0.4G $mcscanshell &";
   system "$mcscan -a $opt{project} > log 2> log2";
}

#system ("perl $mc2 $opt{project}");
#system ("perl $dot $opt{project}");
