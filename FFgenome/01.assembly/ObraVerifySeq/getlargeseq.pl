##get sequence larger than 100kb

$/="\>";
open IN, "VerifySeq20091216.fas" or die "can not open my file";
open OUT, ">largeVerifySeq.fas" or die "can not open my outfile";
while (<IN>){
    chomp $_;
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $head=shift @unit;
    my @array=split(" ",$head);
    $head=$array[0];
    #print "$head\n";
    my $seq=join("",@unit);
    if (length $seq > 100000){
         print OUT ">$head\n$seq\n";
    }
}
close IN;
close OUT;
$/="\n";

