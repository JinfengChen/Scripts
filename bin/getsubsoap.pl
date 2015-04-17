#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"list:s","soap:s","out:s","help");


my $help=<<USAGE;
perl $0 --soap soaprefchr.lst --list revision.mapping.pos --out ../input/soaprefdraw 
--list: list and position of region to fetch soap result
chr11_903308_1634529	chr11	903308	1634529	SD1
chr12_1257880_1838095	chr12	1257880	1838095	SD2
OB_chr04_6039164_6768670	chr04	6039164	6768670	H7
OB_chr04_7027922_7396967	chr04	7027922	7396967	H8
--soap: soap list containg the soap files
--out: output dir of soap
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
}

####read position list
my $reflist=readlist($opt{list});

####read soap list
#/home/jfchen/159/FFproject/FFgenome/01.assembly/PEscafold/input/091023_I77_FC42RRLAAXX_L2_ORYcogDAADLAAPE.soap
my $refsoaplst=readsoaplst($opt{soap});
foreach my $soap (keys %$refsoaplst){
   print "$soap\n";
   my @unit=split("\/",$soap);
   my $out=$opt{out}."/$unit[$#unit]";
   print "$out\n";
   subsoap($soap,$reflist,$out); 
}

#########
#FC42RRLAAXX:2:1:11:1131#0/1	CCCCCCATACTCCAAGAATGTGATTCACGAGGGCCGATAGATCA	BBBBBBBBBBBBBB]U_]^^_[_Q]LTGUV^aa[Raa^aaaa_`	1	a	44	-	chr08	3523550	2A->2C2	T->3C2	44M	2AT40
sub subsoap
{
my ($file,$list,$outfile)=@_;
open IN, "$file" or die "$!";
open OUT, ">$outfile" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my $line=$_;
    my @unit=split("\t",$line);
    my $pair=<IN>;
    chomp $pair;
    my @unit1=split("\t",$pair);
    if (exists $list->{$unit[7]}){
       foreach my $p (@{$list->{$unit[7]}}){
           if ($unit[8] > $p->[1] and $unit[8] < $p->[2] and $unit1[8] > $p->[1] and $unit1[8] < $p->[2]){
               $unit[7]=$p->[0];
               $unit1[7]=$p->[0];
               $unit[8]=$unit[8]-$p->[1]+1;
               $unit1[8]=$unit1[8]-$p->[1]+1;
               my $temp=join("\t",@unit);
               my $temp1=join("\t",@unit1);
               #print "$unit[7]\t$p->[1]\t$p->[2]\n";
               #print "$line\n$pair\n";
               print OUT "$temp\n$temp1\n";
           }
       }
    }
}
close OUT;
close IN;
}


sub readlist
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    push (@{$hash{$unit[1]}},[$unit[0],$unit[2],$unit[3]]);
}
close IN;
return \%hash;
}


sub readsoaplst
{
my ($file)=@_;
my %hash;
open IN, "$file" or die "$!";
while(<IN>){
    chomp $_;
    next if ($_=~/^$/);
    my @unit=split("\t",$_);
    $hash{$unit[0]}=1;
}
close IN;
return \%hash;
}
 
