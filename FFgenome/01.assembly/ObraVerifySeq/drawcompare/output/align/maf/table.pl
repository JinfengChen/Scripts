while (glob("*.maf")){
   push (@file,$_);
}
foreach (@file){
   `perl ../../../bin/statAlign.pl -a $_ -f maf > log`; 
}
