#!/usr/bin/perl
use strict;

##draw the distribution figure for hcluster score

my $file = shift;

my $distribute_svg = "/home/jfchen/FFproject/tools/draw/distribute_svg/distribute_svg.pl";
my $svg2xxx = "/home/jfchen/FFproject/tools/draw/svg2xxx_release/svg2xxx";
my $distribute_pre = "/home/jfchen/FFproject/tools/draw/distribute_assist/bin/distribute_pre.pl";

`more $file | awk '\$1!=\$2' > $file.score`;
`awk '{print \$3}' $file.score  | $distribute_pre --frequency --minborder 0 --maxborder 100 --binsize 1  --header rect --color blue --mark frequency  -Ystart 0 -Yend 10 -Ystep 2  -X \"Hcluster score\" -Y \"Percent in all\" > $file.score.lst`;
`$distribute_svg $file.score.lst $file.score.svg`;
`$svg2xxx $file.score.svg`;
`rm $file.score`;
