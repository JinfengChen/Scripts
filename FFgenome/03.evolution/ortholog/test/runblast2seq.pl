opendir DIR, "./" or die "can not open my dir";
foreach my $file (readdir DIR){
    if ($file=~/(.*)\.fas/){
       print "$file\n";
       push (@seq,$1);
	   push (@name,$1);
       system "formatdb -i $file -p T";
    }
}
close DIR;

for(my $i=0;$i<@seq;$i++){
   for(my $j=0;$j<@seq;$j++){
     system "blastall -p blastp -i $seq[$i].fas -d $seq[$j].fas -o $seq[$i]VS$seq[$j].blast -e 1e-20 -m 8";
     print "$i\t$j\n"; 
   }
}








