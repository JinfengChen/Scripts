#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"outdir:s","help");


my $help=<<USAGE;
perl $0 --outdir chrmcscan_final8_3_last
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

`mkdir $opt{outdir}`;
my $osbed="../input/Os.bed.chr";
my $obbed="../input/Ob.bed.chr";
my $sbbed="../input/Sb.bed.chr";
my $ogbed="../input/Og.bed.chr";
my $oibed="../input/Oi.bed.chr";
my $blast="../input/oryza.blast";

runchr("01","03");
runchr("02","04");
runchr("03","01");
runchr("04","06");
runchr("05","09");
runchr("06","10");
runchr("07","02");
runchr("08","07");
runchr("09","02");
runchr("10","01");
runchr("11","05");
runchr("12","08");

sub runchr
{
my ($oryzachr,$sbchr)=@_;
`cat $osbed/Os$oryzachr.bed $obbed/Ob$oryzachr.bed $oibed/Oi$oryzachr.bed $ogbed/Og$oryzachr.bed $sbbed/Sb$sbchr.bed > chr$oryzachr.bed`;
`cp $blast chr$oryzachr.blast`;
`perl runMcscan.pl -mcl --project chr$oryzachr > chr$oryzachr.log 2> chr$oryzachr.log2`;
#`perl sumBlockV3.pl --block chr$oryzachr.blocks`;
`mv chr$oryzachr.* $opt{outdir}`;
}

