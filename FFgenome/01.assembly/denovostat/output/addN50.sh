n50=`perl ../bin/count_N50.pl $1`
n90=`perl ../bin/count_N90.pl $1`
#n50=`perl ../bin/count_N50.pl ../data/chr04.scafSeq`  
#n90=`perl ../bin/count_N90.pl ../data/chr04.scafSeq`
echo 
echo 
echo 
echo Color:#FA9B9D
echo Mark: N50 size
echo Type:LineBar
echo LineDash:7
echo $n50:0
echo $n50:100
echo 
echo 
echo 
echo Color:#D2691E
echo Mark: N90 size
echo Type:LineBar 
echo LineDash:7 
echo $n90:0 
echo $n90:100


