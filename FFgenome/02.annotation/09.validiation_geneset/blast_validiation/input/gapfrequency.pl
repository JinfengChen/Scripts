

die "Usage: perl gapfrequency.pl name.fas > gaps\n" if (@ARGV < 1);
$/=">";

open IN, "$ARGV[0]" or die "can not open infile";
while (<IN>){
       chomp $_;
       next if ($_ eq "");
       my @unit=split("\n", $_);
       my $head=shift @unit;
       #print "$head\n";
       my $seq=join ("",@unit);
       #next if (length $seq < 100000);
       $seq=~s/\r//g;
       $seq=~s/\n//g;
       $seq=~s/\s//g;
       if ($seq=~/NN/){
            $counter++;
       }
       while ($seq=~/(N+)/g){
           my $len=length $1;
           #print "$len\n";
           my $pos=pos($seq);
           #print "$pos\n"; 
       }

        
}
close IN;
print "$counter\n";
