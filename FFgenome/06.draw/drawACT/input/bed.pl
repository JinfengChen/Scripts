#chr05   22243158        22243194        SRR042529.6 HWI-EAS412:6:1:3:330 length=36/2    25      -

open IN, "$ARGV[0]" or die "$!";
while(<IN>){
    chomp $_;
    my @unit=split("\t",$_);
    my @array=split(" ",$unit[3]);
    print "$unit[0]\t$unit[1]\t$unit[2]\t$array[0]\t$unit[4]\t$unit[5]\n";

}
close IN;
