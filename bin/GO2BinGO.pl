#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"GO:s","help");


my $help=<<USAGE;
Convert iprscan GO to BinGO format that can be used in BinGO modular.
perl $0 --GO 
USAGE


if ($opt{help} or keys %opt < 1){
    print "$help\n";
    exit();
} 

& convert($opt{GO});

#Ob0037g10040.1  3       GO:0003824; catalytic activity; Molecular Function      GO:0008270; zinc ion binding;
sub convert
{
my ($file)=@_;
print "(species=Oryza glaberrima)(type=Biological Process)(curator=GO)\n";
open IN,"$file" or die "$!";
while(<IN>){
     chomp $_;
     my @unit=split("\t",$_);     
     for (my $i=2;$i<@unit;$i++){
         if ($unit[$i]=~/GO:(\d+)\;/){
            print "$unit[0] \= $1\n";
         }
     }
}
close IN;
}
