
my %hash;
open IN, "$ARGV[0]" or die "$!";
while(<IN>){
   chomp $_;
   next unless ($_=~/OBR/);
   my @unit=split("\t",$_);
   my $read=$unit[2]+$unit[3];
   $hash{$unit[0]}=$read;
}
close IN;

my $num=keys %hash;
print "$num\n";
