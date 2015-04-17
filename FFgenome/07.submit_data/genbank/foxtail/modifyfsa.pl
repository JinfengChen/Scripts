## usage: perl modifyfsa.pl -i infile.fas -l locusname -p protein

use Getopt::Long;
my ($infile,$locus,$protein);
GetOptions(
  "infile:s"=>\$infile,
  "locus:s" =>\$locus,
  "protein:s"=>\$protein
);
#print "$infile\n";
die "Usage: perl modifyfsa.pl -i infile -l locusname -p protein\n" unless (-f $infile);
$/=">";
open IN, "$infile" or die "can not open my file";
    while (<IN>){
        $_=~s/\r//g;
        next if (length $_ < 2); 
        my @unit=split("\n",$_);
        my $head=shift @unit;
        my $seq =join("",@unit);
        #print "$head\n";
        $seq=~s/\>//g;
        #print "$seq\n";
        my $species1="Setaria viridis";
        my $species2="Setaria italica";
        my $identify=substr($head,0,1);
        #print "$identify\n";
        if ($identify eq "q"){
            print ">$locus\_$head $species1, $protein ($locus), code: $head [organism=Setaria viridis] [strain=$head]\n$seq\n";
        }elsif($identify eq "n" or $identify eq "g"){
            print ">$locus\_$head $species2, $protein ($locus), code: $head [organism=Setaria italica] [strain=$head]\n$seq\n";
        }else{
            print ">$locus\_$head $species2, $protein ($locus), code: Yugu1 [organism=Setaria italica] [strain=Yugu1]\n$seq\n";
        }
    }
close IN;

