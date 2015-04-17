system "rm *.phr *.psq *pin";
while(glob("*.fas")){
     if ($_=~/(.*)\.fas/){
        push (@species,$1); 
     }
}
my %length;
for(my $i=0;$i<@species;$i++){
   open IN, "$species[$i]\.length" or die "length erros";
        while(<IN>){
            my @unit=split("\t",$_);
            #print "$unit[0]\t$unit[1]\n";
            $length{$unit[0]}=$unit[1];  
        } 
   close IN;
}
for(my $i=0;$i<@species-1;$i++){
   for(my $j=$i+1;$j<@species;$j++){
      my %pair1;
      my %info1;
      my %pair2;
      my %info2;
      my %rate2;
      my %rate1;
      print "$species[$i]\t$species[$j]\n";
      #my $infile1="$species[$i]VS$species[$j]\.blast";
      open IN1, "$species[$i]VS$species[$j]\.blast" or die "can not open my infile1";
           while(<IN1>){
                 #print "$_"; 
                 my @unit=split("\t",$_);
                 unless(exists $pair1{$unit[0]}){
                     my $lengthrate=$unit[3]/$length{$unit[0]};
                     if ($unit[2] >= 30 and $lengthrate >= 0.3){
                         $pair1{$unit[0]}=$unit[1];
                         $info1{$unit[0]}="$unit[2]";
                         $rate1{$unit[0]}=$lengthrate;
                     }         
                 }
           }
      close IN1; 
      open IN1, "$species[$j]VS$species[$i]\.blast" or die "can not open my infile2";
           while(<IN1>){
                 #print "$_"; 
                 my @unit=split("\t",$_);
                 unless(exists $pair2{$unit[0]}){
                     my $lengthrate=$unit[3]/$length{$unit[0]};
                     if ($unit[2] >= 30 and $lengthrate >= 0.3){
                         $pair2{$unit[0]}=$unit[1];
                         $info2{$unit[0]}="$unit[2]";
                     }         
                 }
           }
      close IN1; 
      open OUT, ">$species[$i]2$species[$j].pair" or die "can not open outfile";
           foreach(sort keys %pair1){
                if ($_ eq $pair2{$pair1{$_}}){
                      print OUT "$_\t$pair1{$_}\t$info1{$_}\t$rate1{$_}\n";    
                }
           }
      close OUT;
   }
}
system "find *.pair | xargs wc -l"
