
my @file=glob("chr*");
chomp @file;
foreach(@file){
system("/home/jfchen/FFproject/tools/BEDTools/bin/mergeBed -i $_ > $_.BED.merge");

}
