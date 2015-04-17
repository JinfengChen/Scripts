#!/usr/bin/perl

####Convert fasta alignment in current directory into paml(*.nuc) format
####Run: perl fasta2paml.pl

my $directory=$ARGV[0];
my $file1;
$/=">";

opendir DIR, $directory or die "can not open dir";
          foreach $file1 (readdir DIR) {
			    
                if ($file1=~/(.*)\.fas/) {
                    my $outfile;
                    my $counter;
                    my $length;
                    my @records="";
                    $outfile="$1";
	            open INFILE, "$directory/$file1" or die "can not open my infile";
                    open OUTFILE, ">$directory/$outfile.nuc" or die "can not open my outfile";
		    #open OUTFILE1, ">$directory/$outfile" or die "can not open my outfile1";
                    while (my $line=<INFILE>) {
	                      ++$counter;
	                      chomp $line;
	                      my @word=split("\n",$line); 
	                      my $head=shift @word;
	                      my $seq=join("",@word);
                              $head=~s/\r//g;
                              $seq=~s/\r//g;
                              $seq=~s/\s//g;
                              $seq=~s/\n//g;
			      $seq=~tr/atcg/ATCG/;
	                      $length=length $seq;
	                      #print "$head\n";
			      #print OUTFILE1 "$head\n";
	                      my $record="$head\n$seq\n";
	                      push(@records,$record);
                    }
		    my $number=$counter-1;
                    print OUTFILE "$number   $length\n";
                    foreach my $r(@records) {
	                        if ($r=~/^\s$/) {
		                          next;
	                        }else{
                                  print OUTFILE "$r";
                            }
                    }
                    close INFILE;
                    close OUTFILE;
		    #close OUTFILE1;
                }
                system ("rm $directory/$file1");
          }
closedir DIR;









