for file in *.maf;
do `perl ../../../bin/statAlign.pl -a $file -f maf > log &`;
done;
