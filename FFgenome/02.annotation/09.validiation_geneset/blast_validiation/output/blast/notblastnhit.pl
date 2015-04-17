

my %hit;
open IN, "blastnohit.pep.blast.tab" or die "$!";
while (<IN>){
     next if ($_ eq "" or $_ =~/Query/);
     my @unit=split("\t",$_);
     $hit{$unit[0]}=1;
} 
close IN;


open IN, "RAPnohit" or die "$!";
while (<IN>){
     next if ($_ eq "");
     my @unit=split("\t",$_);
     unless (exists $hit{$unit[0]}){
         print "$_";
     }
}
close IN;
