

my %hash;
open IN, "anchoredscaf.txt" or die "can not open my anchored file";
while (<IN>){
     chomp $_;
     next if ($_ eq "");
     my @unit=split("\t",$_);
     $hash{$unit[0]}=$unit[1];
}
close IN;
open IN, "scaffold" or die "can not open my scaffold";

my $len;
while (<IN>){
     chomp $_;
     next if ($_ eq ""); 
     my @unit=split("\t",$_);
     unless (exists $hash{$unit[0]}){
           $len+=$unit[1];
           print "$_\n"; 
     }
}
print "Total Unanchored Length: $len\n";
close IN;
