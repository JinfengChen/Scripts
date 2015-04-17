
my $directory=$ARGV[0];
open SUMMARY, ">>yn00.summary" or die "can not open my summary";
print SUMMARY "	Nei & Gojobori 1986				Yang & Nielsen 2000				Li W-H 1993 and Pamilo P, Bianchi NO 1993\n";
print SUMMARY "Gene\tKa\tKs\tOmega\tMyr\tKa\tKs\tOmega\tMyr\tKa\tKs\tOmega\tMyr\tAverage Ks\tAverage Year\n";
opendir DIR, $directory or die "can not open dir";
foreach my $file1 (readdir DIR) {	    
                if ($file1=~/(.*)\.nuc/) {
		   my $name=$1;
		   print SUMMARY "$name\t";
		   my $ks;
                   my $ctl="$directory/yn00.$name";
                   open FILE, "/share/raid12/chenjinfeng/FFgenome/evolution/kaks_pairwise/bin/script_yn00/yn00.ctl" or die "can not open my control file";
		   open OUT, ">$ctl" or die "can not open my outfile";
                      while (<FILE>) {
		         if ($_=~/seqfile =/) {
                             print OUT "      seqfile = $file1 * sequence data file name\n"; 
			 }elsif($_=~/outfile =/){
			     print OUT "      outfile = $name.out   * main result file\n";
			 }else{
			     print OUT "$_";
		         }
                      } 
                 close FILE; 
                 close OUT;
	         system "yn00 $ctl";
                 my $out="$directory/$name.out";
                 open RESULT, "$out" or die "can not open my result file";		 
		          while (<RESULT>) {
                          if ($_=~/\(A\) Nei-Gojobori \(1986\) method/) {
			      <RESULT>;
                              <RESULT>;
                              <RESULT>;
                              <RESULT>;
                              <RESULT>;
                              <RESULT>;
			      my $ng=<RESULT>;
			      if ($ng=~/\s+(\d+\.\d+)\((\d+\.\d+) (\d+\.\d+)\)/) {
			          my $year=$3*1000/13;
				  $ks=$ks+$3;
				  print "$name!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";				  
                                  print SUMMARY "$2\t$3\t$1\t$year\t";
			      }
			   }elsif($_=~/\(B\) Yang \& Nielsen \(2000\) method/){
			      <RESULT>;
                              <RESULT>;
                              <RESULT>;
                              <RESULT>;
                              <RESULT>;
                              <RESULT>;
                              <RESULT>;
			      my $yn=<RESULT>;
			      my @unit=split(" ",$yn);
			      my $year=$unit[10]*1000/13;
			      $ks=$ks+$unit[10];
                              print "YN: $yn\t$unit[6]\n";
			      print SUMMARY "$unit[7]\t$unit[10]\t$unit[6]\t$year\t";
			   }elsif($_=~/\(C\) LWL85\, LPB93 \& LWLm methods/){
			      <RESULT>;
                              <RESULT>;
                              <RESULT>;
                              <RESULT>;
                              <RESULT>;
                              <RESULT>;
                              <RESULT>;
			      <RESULT>;
                              <RESULT>;
                              <RESULT>;
                              <RESULT>;
                              <RESULT>;
                              <RESULT>;
                              <RESULT>;
                              <RESULT>;
			      my $lpb=<RESULT>;
			      my @unit=split(" ",$lpb);
                              my $year=$unit[3]*1000/13;
			      $ks=$ks+$unit[3];
			      print SUMMARY "$unit[6]\t$unit[3]\t$unit[9]\t$year\t";
		              my $averageks=$ks/3;
		              my $averageyear=$averageks*1000/13;
		              print SUMMARY "$averageks\t$averageyear\n";
                           }
	                   }
	         close RESULT;
                 system ("rm $directory/$file1 $out $ctl");
                 }
}
closedir DIR;
close SUMMARY;
`rm 2YN.*`;
`rm rst rst1 rub`;
