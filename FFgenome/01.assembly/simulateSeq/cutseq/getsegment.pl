#!/usr/bin/perl
# perl getsegment.pl 1000 20000 chr04.txt

my $start=@ARGV[0];
my $end=@ARGV[1];
my $file=$ARGV[2];
$/=">";			    
my $outfile="$file"."$start:$end"; 
open INFILE, "$file" or die "can not open my infile";
open OUTFILE, ">$outfile.txt" or die "can not open my outfile";
		

           while (my $line=<INFILE>) {
                          if ($line=~/\w+/){
	                      $counter++;
	                      chomp $line;
	                      my @word=split("\n",$line); 
	                      my $head=shift @word;
	                      my $seq=join("",@word);
                              $seq=~tr/atcg/ATCG/;
	                      my $length=$end-$start+1;
			      my $segment=substr($seq,$start,$length);
	                      print OUTFILE ">$outfile\n$segment\n";
						  
                          }
		   }

close INFILE;
close OUTFILE;









