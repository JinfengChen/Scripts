

die "Usage: perl $0 name.fas > gaps\n" if (@ARGV < 1);
$/=">";
open IN, "$ARGV[0]" or die "can not open infile";
while (<IN>){
       chomp $_;
       next if ($_ eq "");
       my @unit=split("\n", $_);
       my $temp=shift @unit;
       my @temp=split(" ",$temp);
       my $head=$temp[0];
       my $seq=join ("",@unit);
       my $seqlen=length $seq;
       $seq=~s/\r//g;
       $seq=~s/\n//g;
       $seq=~s/\s//g;
       while ($seq=~/(N+)/g){
           $counter++;
           my $len=length $1;
           $length+=$len;
       }
}
close IN;
print "Gap Number:$counter\nGap Size:$length\n";
