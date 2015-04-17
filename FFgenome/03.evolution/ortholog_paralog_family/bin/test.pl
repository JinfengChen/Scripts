
open IN, "$ARGV[0]" or die "$!";
while(<IN>){
   chomp $_;
   my $line=$_;
   if ($line=~/OBa_g30$/){
     print "$line\n";
   }
}
close IN;
