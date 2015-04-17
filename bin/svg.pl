#!/usr/bin/perl -w
use strict;
use SVG;
# create an SVG object
my $svg= SVG->new(width=>200,height=>200);    

   my $tag=$svg->text(
          x=> 100, y=> 100,
          style =>{
               'font-size' => 10,
               'text-anchor' => 'start',
               'transform' => 'rotate(-90)'
          }
          #transform => 'rotate(-90)'
       )->cdata("TEST1111111111111");


drawsvg($svg,"testsvg");


 sub drawsvg
 {
 my ($svg,$name)=@_;
 open OUT,">$name\.svg";
 print OUT $svg->xmlify;
 close OUT;
 system("/rhome/cjinfeng/software/tools/draw/svg2xxx_release/svg2xxx -t pdf -m 3024 $name.svg");
 } 

