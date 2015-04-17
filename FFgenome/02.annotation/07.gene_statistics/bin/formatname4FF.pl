
$/=">";
open IN, "$ARGV[0]" or die "$!";
open OUT, ">$ARGV[0].shorthead" or die "$!";
while(<IN>){
    chomp $_;
    next if (length $_ < 2);
    my @unit=split("\n",$_);
    my $head=shift @unit;
    if ($head=~/(\w+) /){
        $head=$1;
    }
    my $seq =join("\n",@unit);
    print OUT ">$head\n$seq\n";
}
close OUT;
close IN;
$/="\n";
