#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"IR:s","help");


my $help=<<USAGE;
perl $0 --IR IRlist.txt

USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

& summary ($opt{IR});


#2  Gene      ::::,::::       147167:DNA/MuDR:2145134:2145304:171,147155:DNA/MuDR:2157500:2157710:209
sub summary
{
my ($file)=@_;
open IN, "$file" or die "$!";
my @sum;
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    if ($unit[3]=~/\w+\,\w+/ and $unit[2]=~/\w+\,\w+/){
       $sum[0]++;
    }elsif($unit[2]=~/\w+\,\w+/){
       $sum[1]++;
    }elsif($unit[3]=~/\w+\,\w+/){
       $sum[2]++;
    }else{
       $sum[3]++;
    }
}
close IN;
print "ALL:$sum[0]\tRef:$sum[1]\tQry:$sum[2]\tNone:$sum[3]\n";
} 
