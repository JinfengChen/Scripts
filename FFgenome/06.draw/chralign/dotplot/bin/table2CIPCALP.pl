#!/usr/bin/perl
use Getopt::Long;

GetOptions (\%opt,"help");


my $help=<<USAGE;
This scripts do the same thing as filter_blast.py. 1. kick self match. 2. Remove multi HSP. But we use CIP/CALP in step 2.
Convert blasttable to Mcscan familiar format filtered by e-value=1e-40 and CIP/CALP. blasttable was generage by blastpaser.pl of BGI.
CIP: Cumulative Identity Percentage. CIP=cumulative percent of sequence identity for all HSP/cumulative aligned length(AL)
CALP: Cumulative Alignment Length Percentage. CALP=AL/query length.

Run: perl table2CIPCALP.pl EC-SC.blasttable > EC-SC.blast

USAGE


if ($opt{help} or @ARGV < 1){
   print "$help\n";
   exit();
}

my $cipcut=0.85;
my $calpcut=0.75;

my %record;
open IN, "$ARGV[0]" or die "$!";
<IN>;
while(<IN>){
   my @unit=split("\t",$_);
   next if ($unit[13] > 1e-40);
   next if ($unit[0] eq $unit[4]); ### kick self match
   my $pair=$unit[0]."-".$unit[4];
   if (exists $record{$pair}){
        my $refarray=$record{$pair};
	push (@$refarray,[@unit]);
	$record{$pair}=$refarray;
   }else{  
	my @array;
	push (@array,[@unit]);
        $record{$pair}=\@array;
   }
}
close IN;

my %inparanoid;
foreach (keys %record){
      #print "$_\n";
      #print Dumper($record{$_}),"\n";
      my $refarray=$record{$_};
      my $num=@$refarray;
      #print "$num\n";
      
          my $queryid =$$refarray[0][0];
	  my $hitid   =$$refarray[0][4];
	  my $querylen=$$refarray[0][1];
	  my $hitlen  =$$refarray[0][5];
          my $alignlen=$$refarray[0][11];
          my $identity=$$refarray[0][8]*$$refarray[0][11];
          my $evalue  =$$refarray[0][13];
      if ($num == 1){
	  my $CIP =$identity/$alignlen;
          my $CALP=$alignlen/$querylen;
          if ($CIP >= $cipcut and $CALP >= $calpcut){
             #$inparanoid{$_}="$queryid\t$hitid\t$querylen\t$hitlen\t$CIP\t$CALP\t$evalue\n";
             $inparanoid{$_}="$queryid\t$hitid\t$evalue\n";
          }
      }else{  
	  my $hspendq  =$$refarray[0][3];
	  my $hspendh  =$$refarray[0][7];
          for(my $i=1;$i<=$num-1;$i++){
                if ($$refarray[$i][2] < $hspendq or $$refarray[$i][6] < $hspendh){
	              next;
		}else{
		      $hspendq=$$refarray[$i][3];
		      $hspendh=$$refarray[$i][7];
                      $alignlen+=$$refarray[$i][11];                
                      $identity+=$$refarray[$i][8]*$$refarray[$i][11];
                      if ($$refarray[$i][13] < $evalue){
                         $evalue=$$refarray[$i][13];
                      }
		}
      
          }
          my $CIP =$identity/$alignlen; 
          my $CALP=$alignlen/$querylen;
          #$inparanoid{$_}="$queryid\t$hitid\t$querylen\t$hitlen\t$CIP\t$CALP\n";
          if ($CIP >= $cipcut and $CALP >= $calpcut){
             $inparanoid{$_}="$queryid\t$hitid\t$evalue\n";
          }
      }
}

foreach (sort keys %inparanoid){
    print "$inparanoid{$_}";
}
