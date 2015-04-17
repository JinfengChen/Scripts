         my %read;
         open IN, "$ARGV[0]" or die "$!";
             while (<IN>){
                  if ($_=~/^@(SRR042638.*)$/){
                       my @unit=split(" ",$1);
                       my $head=$unit[0];
                       if (exists $read{$head}){
                          #print "Dupli\t$head\n";
                          $read{$head}=1;
                       }else{
                          $read{$head}=1;
                          print "$head\n";
                       }
                  }     
             }
         close IN;


