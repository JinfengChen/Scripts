#!/usr/bin/perl
use Data::Dumper
### Convert blasttable to inparanoid format. blasttable was generage by blastpaser.pl of BGI.
### perl table2inparanoid.pl EC-SC.blasttable > EC-SC

my %record;
open IN, "$ARGV[0]" or die "$!";
<IN>;
while(<IN>){
   my @unit=split("\t",$_);
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
	  my $bitscore=$$refarray[0][12];
	  my $querylen=$$refarray[0][1];
	  my $hitlen  =$$refarray[0][5];
          my $longestq=$$refarray[0][3]-$$refarray[0][2]+1;
	  my $longesth=$$refarray[0][7]-$$refarray[0][6]+1;
	  my $matchq  =$$refarray[0][3]-$$refarray[0][2]+1;
	  my $matchh  =$$refarray[0][7]-$$refarray[0][6]+1;
	  my $position="q:$$refarray[0][2]-$$refarray[0][3] h:$$refarray[0][6]-$$refarray[0][7]";
      if ($num == 1){
	  $inparanoid{$_}="$queryid\t$hitid\t$bitscore\t$querylen\t$hitlen\t$longestq\t$longesth\t$matchq\t$matchh\t$position";
      }else{  
          my $hspstartq=$$refarray[0][2];
	  my $hspstarth=$$refarray[0][6];
	  my $hspendq  =$$refarray[0][3];
	  my $hspendh  =$$refarray[0][7];
          for(my $i=1;$i<=$num-1;$i++){
                if ($$refarray[$i][2] < $hspendq or $$refarray[$i][6] < $hspendh){
	              next;
		}else{
		      $hspendq=$$refarray[$i][3];
		      $hspendh=$$refarray[$i][7];
		      $bitscore+=$$refarray[$i][12];
		      $matchq  +=$$refarray[$i][3]-$$refarray[$i][2]+1;
		      $matchh  +=$$refarray[$i][7]-$$refarray[$i][6]+1;
		      $position.="\tq:$$refarray[$i][2]-$$refarray[$i][3] h:$$refarray[$i][6]-$$refarray[$i][7]";
		}
      
          } 
	  $longestq=$hspendq-$hspstartq+1;
	  $longesth=$hspendh-$hspstarth+1;
	  $inparanoid{$_}="$queryid\t$hitid\t$bitscore\t$querylen\t$hitlen\t$longestq\t$longesth\t$matchq\t$matchh\t$position";
      }
}

foreach (sort keys %inparanoid){
    print "$inparanoid{$_}\n";
}
