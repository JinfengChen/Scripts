my %read;
open IN, "$ARGV[0]" or die "$!";
while (<IN>){
   chomp $_;
   my $line=$_;
   my @unit=split("\t",$_);
   if ($unit[0] eq "chr7" and $unit[3]=~/(.*)\/1/){
       if (exists $read{$1}){
            $read{$1}.="\n$unit[0]";
            #print "Dupli\t$1\t$read{$1}\n";
       }else{
            $read{$1}=$unit[0];
            print "$1\n";
       }
   }
}
close IN;

