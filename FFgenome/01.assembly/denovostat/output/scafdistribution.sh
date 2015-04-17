myfile=$1
fastaDeal.pl $myfile -attr len | awk '$1 > 1000' | perl /share/raid12/chenjinfeng/tools/draw/distribute_assist/bin/distribute_pre.pl -frequency -binsize 200 -color red -mark ScafLength -header line -ystart 0 -yend 100 -xstart 0 -xend 4000000 -xscale "/1000"  -note "Scaffold Length Distribution" -x "Scaf Length,Kb" -y "frequency" -fontsize 25 > scaf.lst
./addN50.sh $myfile > addn50.lst
cat scaf.lst addn50.lst > scaflength.lst
rm scaf.lst addn50.lst
/share/raid12/chenjinfeng/tools/draw/distribute_assist/bin/distribute_svg.pl scaflength.lst scaflength.svg
