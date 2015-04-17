#!/usr/bin/perl


my $file1=$ARGV[0];
$/=">";
if ($file1=~/\w+/) {
                    my $outfile;
                    my $counter=1;
                    my $length;
                    my @records="";
		    my @name;
		    my %sequence;
                    $outfile="$file1.phy";
		    open INFILE, "$file1" or die "can not open my infile";
                    open OUTFILE, ">$outfile" or die "can not open my outfile";
                    while (my $line=<INFILE>) {
	                      chomp $line;
	                      my @word=split("\n",$line); 
	                      my $head=shift @word;
                              $head=$1 if ($head=~/(.*?)\s+/);
                              my $seq=join("",@word);
                              $head=~s/\r//g;
                              $seq=~s/\r//g;
                              $seq=~s/\s//g;
                              $seq=~s/\n//g;
                              $length=length $seq;
                              my $record="$counter  $seq\n";
                              if ($head=~/\w+/) {
                                if ($head=~/(.*)\/(.*)/){
                                    $head=$1;
                                }
                                if (length $head > 50) {
                                    $head=substr($head,-20);
                                }
                                push (@name,$head);
                                $sequence{$head}=$seq;
                                push(@records,$record);
                                ++$counter;
                              }###if 
                    }
		    my $number=$counter-1;
                    print OUTFILE "$number   $length\n";
                    foreach my $r(@name) {
	                        if ($r=~/^\s$/) {
		                  next;
	                        }else{
			          printf OUTFILE "%-10s","$r";
                                  print OUTFILE "  $sequence{$r}\n";
                                }
                    }
                    close INFILE;
                    close OUTFILE;
}









