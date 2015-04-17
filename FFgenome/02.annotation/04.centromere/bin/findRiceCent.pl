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


foreach (sort keys %scaf2cent){
       my @array= sort {$a <=> $b} @{$scaf2cent{$_}};
       for (my $i=0;$i<2;$i++){
          shift @array; ### delete the fisrt ten;
          #pop @array;   ### delete the last ten; 
       }
       my $a=$array[0];
       my $b=$array[$#array];
       #my $start=join(",", @{$scaf2cent{$_}});
       my $start=join(",", @array);
       if (@array >=3 ){
          print "$_\t$a\t$b\t$start\n";
       }
}


