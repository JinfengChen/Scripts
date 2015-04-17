### find centromere scaffold if it have more than 10 copies of cent repeat that > 100bp

## Usage: perl findcentscaf.pl ../data/cent2scafblast2 > centromere.txt &
use FindBin qw ($Bin);
my $input="$Bin/../data";
my %scaf2cent;
open IN, "$ARGV[0]" or die "$!";
while (<IN>){
   chomp $_;
   my @unit=split("\t",$_);
   if ($unit[2] >= 70 and $unit[3] >= 100){
      if (exists $scaf2cent{$unit[1]}){
           my $temp=$scaf2cent{$unit[1]};
           push (@$temp,$unit[8]);
           $scaf2cent{$unit[1]}=$temp;
      }else{
           my @temp;
           push (@temp,$unit[8]); 
           $scaf2cent{$unit[1]}=\@temp;
      }
   }   

}
close IN;

my %scaflen;
open IN, "$input/scaffold" or die "$!";
while (<IN>){
   chomp $_;
   next if ($_ eq "");
   my @unit=split("\t",$_);
   $scaflen{$unit[0]}=$unit[1];
}
close IN;

foreach (keys %scaf2cent){
       my $count=@{$scaf2cent{$_}};
       my $start=join(",",@{$scaf2cent{$_}});
       print "$_\t$count\t$scaflen{$_}\t$start\n";
}


