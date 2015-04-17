#!/usr/bin/perl
use Getopt::Long;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);


GetOptions (\%opt,"mping:s", "confident", "project:s","help");


my $help=<<USAGE;
perl $0 --mping --confident
Convert reloceTE confident/non_ref to gff
input:
mping	TTA	HEG4_2	Chr1:588832..588834	+	5	5	GCATCTTCTTGCATTGGTAGCAAGAAAACGGCAACATATGGGCCTCCGATGGAAGCCACGTCCTGTCCAATTCCCAAACTAACCCA
or
mPing   GGC     HEG4    Chr3    10345176..10345178      -       T:3     R:2     L:1

output:
Chr1    HEG4_2.3.mPing.all      Non_reference   588832  588834  .       .       .       ID=mPing_1;
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

$opt{project} ||= "RelocaTE";

if ($opt{confident}){
    readtable($opt{mping});
}else{
    #print 'all non_ref\n';
    readtable2($opt{mping});
}

#mPing   GGC     HEG4    Chr3    10345176..10345178      -       T:3     R:2     L:1
sub readtable2
{
my ($file)=@_;
my %hash;
my $gff=$file;
$gff=~s/\.txt/\.gff/;
#my $prefix=basename($file,".txt");
my $cout;
open OUT, ">$gff";
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    #next if ($_=~/^$/ or $_=~/single/ or $_=~/insuff/);
    next if ($_=~/^TE/);
    my @unit=split("\t",$_);
    $cout++;
    my ($chr,$start,$end);
    $chr = $unit[3];
    if ($unit[4]=~/(\d+)\.\.(\d+)/){
       $start=$1;
       $end  =$2;
    }
    #print "Chr:$chr $start $end\n";
    $unit[7] =~ s/\D+//g;
    $unit[8] =~ s/\D+//g;
    print OUT "$chr\t$unit[2]\t$opt{project}\t$start\t$end\t.\t.\t.\tID=mPing_$cout;TSD=$unit[1];Right=$unit[7];Left=$unit[8];\n";
    #print OUT "$unit[3]\t$unit[2]\t$opt{project}\t$start\t$end\t.\t.\t.\tID=mPing_$cout;\n";
}
close IN;
return \%hash;
}


#mping	TTA	HEG4_2	Chr1:588832..588834	+	5	5	GCATCTTCTTGCATTGGTAGCAAGAAAACGGCAACATATGGGCCTCCGATGGAAGCCACGTCCTGTCCAATTCCCAAACTAACCCA
sub readtable
{
my ($file)=@_;
my %hash;
my $gff;
$gff=~s/\.txt/\.gff/;
#my $prefix=basename($file,".txt");
my $cout;
open OUT, ">$gff";
open IN, "$file" or die "$!";
<IN>;
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $cout++;
    my ($chr,$start,$end);
    if ($unit[3]=~/(\w+)\:(\d+)\.\.(\d+)/){
       $chr  =$1;
       $start=$2;
       $end  =$3;
    }
    print OUT "$chr\t$unit[2]\t$opt{project}\t$start\t$end\t.\t.\t.\tID=mPing_$cout;TSD=$unit[1];Right=$unit[5];Left=$unit[6];\n";
    #print OUT "$unit[3]\t$unit[2]\t$opt{project}\t$start\t$end\t.\t.\t.\tID=mPing_$cout;\n";
}
close IN;
return \%hash;
}
