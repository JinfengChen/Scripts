mkdir maf
mkdir pdf
for file in *.fas
do ./fastaDeal.pl -cuts 1 $file;
echo $file
lastz ./$file.cut/$file.1  ./$file.cut/$file.2 --identity=95..100 --strand=both --gfextend --gapped  --chain --format=maf --rdotplot=./pdf/4r > ./maf/$file.maf;
#perl statAlign.pl -a ./maf/$file.maf -f maf > log 
#perl unalign.pl -a ./maf/$file.maf -f maf >> unalign.table
#cat lastzplot.r | R --vanilla --slave;
#mv ./pdf/4r.pdf ./pdf/$file.pdf;
#rm ./pdf/4r;
done;
