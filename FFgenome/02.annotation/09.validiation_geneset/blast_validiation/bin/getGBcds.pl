#!/usr/bin/perl



$/="//\n";
open IN, "$ARGV[0]" or die "$!";
open OUT, ">$ARGV[0].seq" or die "$!";
while (<IN>){
     my $counter;
     my @unit=split("\n",$_);
     my @array=split(" ",$unit[0]); 
     if ($_=~/(LOCUS\s+.*\nFEATURES\s+Location\/Qualifiers\n)(.*)\nORIGIN(.*)/s){
          my $title=$1;
          my $feature=$2;
          my $sequence=$3;
          $sequence=~s/\d//sg;
          $sequence=~s/\s//sg;
          $sequence=~s/\///sg;
          my $length=length $sequence;
          print "\n",$length,"\n";
          while ($feature=~/CDS\s+(.*?)\//gs){
             $counter++;
             my $name=$array[1]."_"."$counter";
             my $cds=$1;
             my $strand;
             my $codingseq;
             if ($cds=~/complement/){
                $strand=0;
             }else{
                $strand=1;
             }
              
             $cds=~s/[a-zA-Z]//sg;
             $cds=~s/\(|\)|\>|\<|\s//sg;
             print "$cds\n";
             my @exons=split(",",$cds);
             foreach my $exon (@exons){
                 chomp $exon;
                 $exon=~/(\d+)\D+(\d+)/;
                 my $subseq=substr($sequence,$1-1,$2-$1+1);
                 $codingseq.=$subseq;
             } 
             if ($strand==1){
                 print OUT ">$name\n$codingseq\n";
             }else{
                 my $revcom=reverse $codingseq;
                 $revcom=~tr/ATCGatcg/TAGCtagc/;
                 print OUT ">$name\n$revcom\n";
             }
          }
     }else{
          print "Non match\n";
     }
}
close OUT;
close IN;
