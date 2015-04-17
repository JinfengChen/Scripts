#!/usr/bin/perl
use strict;

die "Input file not specified\n" if (@ARGV<2);

chomp $ARGV[0];
chomp $ARGV[1];
print "$ARGV[0]\n$ARGV[1]\n";
`/share/raid12/chenjinfeng/tools/bin/fastaDeal.pl --pa $ARGV[0] ../input/FFversion2/super-scaffold.fa > $ARGV[0]`;
`/share/raid12/chenjinfeng/tools/bin/fastaDeal.pl --pa $ARGV[1] ../input/FFversion2/super-scaffold.fa > $ARGV[1]`;
`bl2seq -i $ARGV[0] -j $ARGV[1] -p blastn -o ../output/$ARGV[0]VS$ARGV[1].blast -e 1e-100 -D 1`;
#`perl hitbes.pl $ARGV[0] $ARGV[1]`;
#`lastz $ARGV[0] $ARGV[1] --strand=plus --gfextend --gapped  --chain --format=axt --rdotplot=4r > 4axt`;
#`cat lastzplot.r | R --vanilla --slave;`;
#`mv 4r.pdf ../output/$ARGV[0]VS$ARGV[1].pdf`;
`rm $ARGV[0] $ARGV[1]`;
#`rm 4r 4axt`;
