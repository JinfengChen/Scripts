

die "Usage: perl gapfrequency.pl name.fas > gaps\n" if (@ARGV < 1);
$/=">";
open OUT, ">gapflanking.fa" or die "$!";
open IN, "$ARGV[0]" or die "can not open infile";
while (<IN>){
       chomp $_;
       next if ($_ eq "");
       my @unit=split("\n", $_);
       my $temp=shift @unit;
       my @temp=split(" ",$temp);
       my $head=$temp[0];
       #print "$head\n";
       my $seq=join ("",@unit);
       my $seqlen=length $seq;
       #next if (length $seq < 100000);
       my $dir="./gaps";
       system ("mkdir $dir") unless (-e $dir);
       open OUT1, ">./gaps/$head\.gaps" or die "$!";
       print OUT1 "$head\t$seqlen\n";
       $seq=~s/\r//g;
       $seq=~s/\n//g;
       $seq=~s/\s//g;
       while ($seq=~/(N+)/g){
           $counter++;
           my $len=length $1;
           print "$len\n";
           my $pos=pos($seq);
           #print "$pos\n";
           my $gaps=$pos-$len;
         if ($len <= 200){
           print OUT1 "$gaps\t$pos\t$len\n";
           my ($start1,$start2);
           if ($pos+200 > length $seq){
              next;
           }else{
              $start2=$pos;
           }
           if ($pos-$len-200 < 0){
              next;
           }else{
              $start1=$pos-$len-200;
           }
           my $id5=$head."_".$pos."_5";
           my $id3=$head."_".$pos."_3";
           my $gap5=substr($seq,$start1,200); 
           my $gap3=substr($seq,$start2,200);
           print OUT ">$id5\n$gap5\n";
           print OUT ">$id3\n$gap3\n";
         }
       }
       close OUT1;
        
}
close IN;
close OUT;
