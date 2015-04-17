open IN, "$ARGV[0]" or die "$!";
my %hash;
while(<IN>){
  next if ($_ eq "");
  $_=~s/\r//g;
  $_=~s/\n//g;
  $hash{$_}=1; 
}
close IN;

foreach (sort keys %hash){
   print  "$_\n";
}
