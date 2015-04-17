

my %hashf;
my %hashfn;
my %hashrn;
my %hashr;
open BES, "../brachyantha.bes.blastm8" or die "can not open my blastout";
while (<BES>){
      my @unit=split("\t",$_);
      if ($unit[0]=~/OB__B(\w+).f/){
         if (exists $hashf{$1}){
            $hashfn{$1}+=1;
            next;
         }else{
            $hashfn{$1}+=1;   
            $hashf{$1}="$unit[1]\t$unit[2]\t$unit[8]\t$unit[9]";
         }
      }
      if ($unit[0]=~/OB__B(\w+).r/){
         if (exists $hashr{$1}){
            $hashrn{$1}+=1;
            next;
         }else{
            $hashrn{$1}+=1;
            $hashr{$1}="$unit[1]\t$unit[2]\t$unit[8]\t$unit[9]";
         }
      } 
}
close BES;

open IN, "sortbycontig" or die "can not open my file 2";
    while (<IN>){
         my @unit=split("\t",$_);
         if (exists $hashf{$unit[0]} and $hashfn{$unit[0]} < 3 and exists $hashr{$unit[0]} and $hashrn{$unit[0]} < 3){
             print "$unit[0]\t$unit[1]\t$hashf{$unit[0]}\t$hashfn{$unit[0]}\t$hashr{$unit[0]}\t$hashrn{$unit[0]}\n"; 
         }
    }
close IN;



