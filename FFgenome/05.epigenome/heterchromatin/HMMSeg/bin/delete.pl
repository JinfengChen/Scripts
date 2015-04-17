
chomp $ARGV[0];
my @file=glob("$ARGV[0]/*local*");
foreach(@file){
   `rm $_`;
}
