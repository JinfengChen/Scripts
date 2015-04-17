#!/usr/bin/perl

# this program is to get a segment of sequence in a certain range that given by user.
# for examle, ./get_subseq.pl 1201 1459 filename newfilename, then, a fasta file will be generated containing the sequence 1201-1459.

my $begin=$ARGV[0]-1;
my $end=$ARGV[1]-1;
my $filename=$ARGV[2];
my $newfile=$ARGV[3];
my $head="";
my $seq="";
my $segment="";

# get the head line and sequence from file through subroutine get_fasta_seq.
($head, $seq)=get_fasta_seq ($filename);

#pick out the segment that user want to find.
$segment=substr ($seq, $begin, $end-$begin+1);

#write the segment to a newfile.
chomp ($head);
my $newhead="$head.segment\n";
open (WRITE, ">$newfile") or die "can not open file";
print WRITE "$newhead$segment\n";
close WRITE;

#print the result file out on screen.
print "$newhead$segment\n";



########subroutines#################################

sub get_fasta_seq{
	
my ($file)=@_;
my $head="";
my $seq="";

open (GET_FASTA_SEQ, "$file") or die "can not open file: $file";

while (my $line=<GET_FASTA_SEQ>) {

      if ($line=~/^\s$/) {
		  next;
	  }elsif ($line=~/^\>/) {
		  $head=$line;
	  }else {
          $seq.=$line;
      }
}
$seq=~s/[\s\n0-9]//g;
close GET_FASTA_SEQ;
return ($head, $seq);
}

