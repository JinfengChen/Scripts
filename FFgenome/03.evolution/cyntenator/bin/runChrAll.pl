#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"help");


my $help=<<USAGE;
perl $0 
USAGE


if ($opt{help}){
    print "$help\n";
    exit();
} 

my $osbed="../input/Os.chr";
my $obbed="../input/Ob.chr";
my $sbbed="../input/Sb.chr";
my $blast="../input/3way.blast";

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
`cp $osbed/Os$oryzachr.txt ./`;
`cp $obbed/Ob$oryzachr.txt ./`;
`cp $sbbed/Sb$sbchr.txt ./`;
`cp $blast chr$oryzachr.blast`;
`./cyntenator -t "((Os$oryzachr.txt Ob$oryzachr.txt) Sb$sbchr.txt)" -h blast chr$oryzachr.blast  > chr$oryzachr.log 2> chr$oryzachr.log2 &`;
}

