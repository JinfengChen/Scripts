while (glob("*.maf")){
   `perl ../../../bin/statAlign.pl -a $_ -f maf > log &`; 
}
