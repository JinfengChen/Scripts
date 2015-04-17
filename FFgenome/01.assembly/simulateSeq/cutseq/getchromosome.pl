#!/usr/bin/perl

my $chromosome=@ARGV[0];
my $file1=@ARGV[1];
$/=">";			    
 
open INFILE, "$file1" or die "can not open my infile";
open OUTFILE, ">$chromosome.txt" or die "can not open my outfile";
		

           while (my $line=<INFILE>) {

	                      $counter++;
	                      chomp $line;
	                      my @word=split("\n",$line); 
	                      my $head=shift @word;
	                      my $seq=join("\n",@word);
						  $seq=~tr/atcg/ATCG/;
	                      $length=length $seq;
	                      if ($head=~/^$chromosome\|/) {
							  print OUTFILE ">$head\n$seq\n";
	                      }
						  

		   }

close INFILE;
close OUTFILE;









