## read .fpc file and calculate contig length for contigs
## usage: perl contiglength.pl *.fpc

if (@ARGV<1){
      die "Usage: perl contiglength.pl *.fpc \n";
      exit(0);
}

my $file=$ARGV[0];
my $contig=422;
for(my $i=1;$i<=$contig;$i++){
   my $ctg="ctg".$i;
   my $chr;
   if ($i<=37){
      $chr=1;
   }elsif ($i<=76){
      $chr=2;
   }elsif ($i<=113){
      $chr=3;
   }elsif ($i<=149){
      $chr=4;
   }elsif ($i<=181){
      $chr=5;
   }elsif ($i<=215){
      $chr=6;
   }elsif ($i<=248){
      $chr=7;
   }elsif ($i<=273){
      $chr=8;
   }elsif ($i<=293){
      $chr=9;
   }elsif ($i<=307){
      $chr=10;
   }elsif ($i<=328){
      $chr=11;
   }elsif ($i<=355) {
      $chr=12;
   }else{
      $chr=0;
   }
    
   my $length=`grep "$ctg" $file | grep "Ends Right" | cut -d " " -f 5 | sort -n | tail -n 1`;
   chomp $length;
   if ($ctg=~/ctg(\d+)/){$ctg=$1}
   $chr="chr".$chr;
   print "$ctg\t$length\t$chr\n";
}
