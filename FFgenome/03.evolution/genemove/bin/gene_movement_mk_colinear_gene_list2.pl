#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"help");


my $help=<<USAGE;
check colinear_gene_list preduced by gene_movement_mk_colinear_gene_list.pl, find non-colinear gene that move within regions.
print "1.arg = global synteny table"
bdid    bdpos   osid    ospos   sbid    sbpos   annotation
USAGE


if ($opt{help} or @ARGV < 1){
    print "$help\n";
    exit();
}


open(IN,"$ARGV[0]");
@global = <IN>;
close(IN);
open (COL,">colinear_genes");
open (NON,">non_colinear_genes");
for($i=0;$i<(@global);$i++) {
   ($osch1,$osgen1,$sbch1,$sbgen1) = ($global[$i] =~ /^\S+\s+\S+\s+Os(\d+)g(\d+)\s+\S+\s+Sb(\d+)g(\d+)\s/);
   $osgen1 =~ s/^[0]+//g;
   $sbgen1 =~ s/^[0]+//g;
   $count_os = 0;
   $count_sb = 0;
   # check low numbers of neighbour genes
   for($j=($i-5);$j<$i;$j++) {
      ($osch2,$osgen2,$sbch2,$sbgen2) = ($global[$j] =~ /^\S+\s+\S+\s+Os(\d+)g(\d+)\s+\S+\s+Sb(\d+)g(\d+)\s/);
      $osgen2 =~ s/^[0]+//g;
      $sbgen2 =~ s/^[0]+//g;
      $diff_os = abs($osgen1-$osgen2);
      $diff_sb = abs($sbgen1-$sbgen2);
      if (($diff_os < 400) || ($osch2 ne $osch1)) {
         $count_os++;
      }
      if (($diff_sb < 400) || ($sbch2 ne $sbch1)) {
         $count_sb++;
      }
      # print "$osch2\t$osgen2\t$sbch2\t$sbgen2\t$count_os\t$count_sb\t$osgen1-$osgen2 = $diff_os\n";
   }
   #print "$global[$i]";
   # check higher numbers of neighbour genes
   for($j=($i+1);$j<$i+5;$j++) {
      ($osch2,$osgen2,$sbch2,$sbgen2) = ($global[$j] =~ /^\S+\s+\S+\s+Os(\d+)g(\d+)\s+\S+\s+Sb(\d+)g(\d+)\s/);
      $osgen2 =~ s/^[0]+//g;
      $sbgen2 =~ s/^[0]+//g;
      $diff_os = abs($osgen1-$osgen2);
      $diff_sb = abs($sbgen1-$sbgen2);
      if (($diff_os < 400) || ($osch2 ne $osch1)) {
         $count_os++;
      }
      if (($diff_sb < 400) || ($sbch2 ne $sbch1)) {
         $count_sb++;
      }
      # print "$osch2\t$osgen2\t$sbch2\t$sbgen2\t$count_os\t$count_sb\t$osgen1-$osgen2 = $diff_os\n";
   }
   if (($count_sb < 4) || ($count_os < 4)) {
      print NON "$global[$i]";
      print "non colinear \t$count_os\t$count_sb \t++++++++++++++++++++++++++++++++\n";
   }else {
      print COL "$global[$i]";
   }
}
# +++++++++++++++++++++++++++++++++++++++++++
sub check_num {
my($k);
} 


