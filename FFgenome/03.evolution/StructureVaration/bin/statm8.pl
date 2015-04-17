#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"m8:s","inf:s","project:s","help");


my $help=<<USAGE;
perl $0 --m8 --inf
--inf FPC infromation of BES and clone.
--m8 blast m8
--project: 
USAGE

if ($opt{help} or keys %opt < 1){
   print "$help\n";
   exit ();
}

my $fpcinf=readfpcinf($opt{inf});
my $title=$1 if ($opt{m8}=~/(.*)\.m8/);
$opt{project}||=$title;
`perl /home/jfchen/FFproject/tools/bin/bestAlign.pl $opt{m8} --cutoff 1e-5 > $opt{project}.best.m8`;


#####stat best hit psl
my $stat=statpsl("$opt{project}.best.m8");
my ($total,$pe,$se);
foreach my $clone (keys %$stat){
    #print "$clone\t$stat->{$clone}->[0]\t$stat->{$clone}->[1]\n";
    $total++;
    if ($stat->{$clone}->[0] == $stat->{$clone}->[1]){
       $pe++;
    }else{
       $se++;
    }
}

####print summary
my $singleBES=$fpcinf->{$opt{project}}->[5]+$fpcinf->{$opt{project}}->[6]-2*$fpcinf->{$opt{project}}->[7];
print "Species\tGenome\tDesinature\tClone with PE\tClone with SE\tClone with PE mapped\tClone with SE mapped\n";
print "$fpcinf->{$opt{project}}->[0]\t$fpcinf->{$opt{project}}->[1]\t$fpcinf->{$opt{project}}->[2]\t$fpcinf->{$opt{project}}->[7]\t$singleBES\t$pe\t$se\n";




########################################
sub readfpcinf
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $name=$1 if ($unit[0]=~/O\. (\w+)/);
    $hash{$name}=[@unit];
}
close IN;
return \%hash;
}

sub statpsl
{
my ($file)=@_;
my %clone;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    my $name=$unit[0];
    if ($name=~/\w+\_\_\w+?(a\d+\D+\d+)\.(\w+)/){
       my $clone=$1;
       my $direction=$2;
       #print "$clone\t$direction\n";
       if ($direction eq "r"){
           $clone{$clone}->[1]=1;
       }elsif($direction eq "f"){
           $clone{$clone}->[0]=1;
       }
    }
}
close IN;
return \%clone;
}

