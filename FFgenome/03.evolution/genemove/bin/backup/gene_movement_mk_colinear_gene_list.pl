#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"help");


my $help=<<USAGE;
print "1.arg = global synteny table"
bdid	bdpos	osid	ospos	sbid	sbpos	annotation	
USAGE


if ($opt{help} or @ARGV < 1){
    print "$help\n";
    exit();
}


open(IN,"$ARGV[0]");
while (<IN>) {
      ($bdid,$bdpos,$osid,$ospos,$sbid,$sbpos,$rem) = /^(\S+)\s+(\S+)\s+(Os\S+)\s+(\S+)\s+(Sb\S+)\s+(\S+)\s+\d+\s+(\S.*)/;
      if ($bdid) {
         ($bdch,$bdgen) = ($bdid =~ /..(\d+)g(\d+)/);
         ($osch,$osgen) = ($osid =~ /..(\d+)g(\d+)/);
         ($sbch,$sbgen) = ($sbid =~ /..(\d+)g(\d+)/);
         $osch =~ s/^0//;
         $sbch =~ s/^0//;
         $bdch =~ s/^0//;
         $bdgen =~ s/^[0]+//;
         $osgen =~ s/^[0]+//;
         $sbgen =~ s/^[0]+//;
         &adjust_chr;
         &match_chr;
         $all_synt = 0;
         #if (($os_sb_synt ==1) && ($os_bd_synt ==0) && ($bd_sb_synt ==0)) {
         if ($os_bd_synt ==1) {
            print "$bdid\t$bdpos\t$osid\t$ospos\t$sbid\t$sbpos\t$rem\n";
         }
     }
}
# adjust ch numbers to compensate for Sbic chromosome fusions
sub adjust_chr {
   # assign virtual Sbic chromosomes 11 + 12
   if (($sbch == 1) && ($sbgen > 16880) && ($sbgen < 31290)) {
      $sbch = 11;
   }
   if (($sbch == 2) && ($sbgen > 11385) && ($sbgen < 32815)) {
      $sbch = 12;
   }
   # assign virtual Bdis chromosomes 6 - 12
   #Bd1
   if (($bdch == 1) && ($bdgen >= 16450) && ($bdgen < 29130)) {
      $bdch = 7;
   }
   if (($bdch == 1) && ($bdgen >= 52490) && ($bdgen <= 59820)) {
      $bdch = 7;
   }
   if (($bdch == 1) && ($bdgen >= 29130) && ($bdgen < 52490)) {
      $bdch = 6;
   }
   #Bd2
   if (($bdch == 2) && ($bdgen >= 14080) && ($bdgen < 40150)) {
      $bdch = 12;
   }
   # Bd3
   if (($bdch == 3) && ($bdgen >= 12460) && ($bdgen < 20730)) {
      $bdch = 8;
   }
   if (($bdch == 3) && ($bdgen >= 34610) && ($bdgen <= 43080)) {
      $bdch = 8;
   }
   if (($bdch == 3) && ($bdgen >= 20730) && ($bdgen < 34610)) {
      $bdch = 10;
   }
   # Bd4
   if (($bdch == 4) && ($bdgen >= 8170) && ($bdgen < 8970)) {
      $bdch = 9;
   }
   if (($bdch == 4) && ($bdgen >= 26920) && ($bdgen < 39020)) {
      $bdch = 9;
   }
   if (($bdch == 4) && ($bdgen >= 8970) && ($bdgen < 26920)) {
   $bdch = 11;
   }
   # define syntenic pairs
   @os_sb_pairs = qw (1-3 2-4 3-1 4-6 5-9 6-10 7-2 8-7 9-12 10-11 11-5 12-8);
   @os_bd_pairs = qw (1-2 2-3 3-1 4-5 5-12 6-6 7-7 8-8 9-9 10-10 11-11 12-4);
   @bd_sb_pairs = qw (1-1 2-3 3-4 5-6 6-10 7-2 8-7 9-12 10-11 11-5 12-9);
} # end adjust_chr -----------------------------------------
# check for chromosome syntheny +++++++++++++++++++++++++++++++
sub match_chr {
    $os_sb_synt = 0;
    SYNT: foreach $pair (@os_sb_pairs) {
             ($os,$sb) = ($pair =~ /(\d+)-(\d+)/);
             if (($osch == $os) && ($sbch == $sb)) {
                $os_sb_synt = 1;
                last SYNT;
             }
          }
    $os_bd_synt = 0;
    SYNT: foreach $pair (@os_bd_pairs) {
             ($os,$bd) = ($pair =~ /(\d+)-(\d+)/);
             if (($osch == $os) && ($bdch == $bd)) {
                $os_bd_synt = 1;
                last SYNT;
             }
          }
          $bd_sb_synt = 0;
    SYNT: foreach $pair (@bd_sb_pairs) {
             ($bd,$sb) = ($pair =~ /(\d+)-(\d+)/);
             if (($bdch == $bd) && ($sbch == $sb)) {
                 $bd_sb_synt = 1;
                 last SYNT;
             }
          }
} # end match_chr ------------------------------------------- 


