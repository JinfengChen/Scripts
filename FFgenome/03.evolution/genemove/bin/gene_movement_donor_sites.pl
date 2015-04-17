#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"help");


my $help=<<USAGE;
print "1.arg = all_grass_colinear_genes_list_final"
print "2.arg = <species>_non_colinear_genes_list"
USAGE


if ($opt{help} or @ARGV < 1){
    print "$help\n";
    exit();
}

open (IN,"$ARGV[1]");
@global = <IN>;
close(IN);
#open(OUT,">donor_site_list");
my %hash;
my %find;
open(IN,"$ARGV[0]");
while (<IN>) {
     ($bdid,$bdpos,$osid,$ospos,$sbid,$sbpos,$rem) = /^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S.*)/;
     if ($bdid) {
         ($osch,$osgen) = ($osid =~ /Os(\d+)g(\d+)/);
         $osgen =~ s/^[0]+//g;
         ($sbch,$sbgen) = ($sbid =~ /Sb(\d+)g(\d+)/);
         $sbgen =~ s/^[0]+//g;

         $donor_yes = 0;
         $entry = "TEST:$bdid\t$bdpos\t$osid\t$sbid\t$rem";
         open (IN2,"$ARGV[1]");
            while (<IN2>) {
                 ($bdid2,$bdpos2,$osid2,$ospos2,$sbid2,$sbpos2,$rem2) = /^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S.*)/;
                 if ($bdid2) {
                    ($osch2,$osgen2) = ($osid2 =~ /Os(\d+)g(\d+)/);
                    $osgen2 =~ s/^[0]+//g;
                    ($sbch2,$sbgen2) = ($sbid2 =~ /Sb(\d+)g(\d+)/);
                    $sbgen2 =~ s/^[0]+//g;
                    $diff_os = abs($osgen-$osgen2);
                    $diff_sb = abs($sbgen-$sbgen2);
                    my $hit = '';
                    if ($osid eq $osid2) { $hit = "+";}  ## donor site has homolog gene in rice
                    if ($sbid eq $sbid2) { $hit .= "-";} ## donor site has homolog gene in sorghum
                    ### all line in $ARGV[1] are colinear in os-sb but moved in brd, so we only need to find the region.
                    ### +- mean donor site have all three gene in colinear, the movement is a copy-paste
                    ### + or - mean donor site may have all three gene in colinear, but not sure
                    ### [blank] mean donor site do not have brd homolog gene, the movement may be a true morement
                    if ((($osch2 eq $osch) && ($diff_os < 1000)) || (($sbch2 eq $sbch) && ($diff_sb < 1000))) {
                       $donor_yes = 1;
                       my $entry1= "$bdid2\t$bdpos2\t$osid2\t$sbid2\n";
                       $hash{$entry1}.="$entry\t$hit\n";
                       $find{$bdid2}.=$hit;
                    }
                }
            }
        close(IN2);
        if ($donor_yes ==1) {
           #print "$entry\n";
        }
     }
} 
my %test;
my $total=@global;
my $find =keys %hash;
foreach my $hit (keys %find){
   if ($find{$hit}=~/\+/ and $find{$hit}=~/\-/){
      $find{1}++;
   }elsif($find{$hit}=~/\+/ or $find{$hit}=~/\-/){
      $find{2}++;
   }else{
      $find{3}++;
   }
}
print "Total movement:$total\nFind Donor:$find\nType1 +-:$find{1}\tType2 +/-:$find{2}\tType3 [blank]:$find{3}\n";
foreach my $g (keys %hash){
   print "$g\n$hash{$g}\n";
}

