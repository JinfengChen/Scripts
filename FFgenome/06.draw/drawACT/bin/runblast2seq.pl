opendir DIR, "./" or die "can not open my dir";
foreach my $file (readdir DIR){
    if ($file=~/(.*)\.fasta/){
       print "$file\n";
       push (@seq,$1);
	   push (@name,$1);
       system "formatdb -i $file -p F";
    }
}
close DIR;

for(my $i=0;$i<@seq;$i++){
   for(my $j=0;$j<@seq;$j++){
     system "blastall -p blastn -U -i $seq[$i].fasta -d $seq[$j].fasta -o $seq[$i]VS$seq[$j].blast -e 1e-10";
     print "$i\t$j\n"; 
   }
}








