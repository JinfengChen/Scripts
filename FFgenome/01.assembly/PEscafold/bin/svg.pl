#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use SVG;


my $svg=SVG->new(width=>800,height=>200);
      my $path=$svg->group(
         style=>{fill=>'none','stroke-wide'=>2}
      );

      my $tag=$svg->line(
         x1=>100,y1=>100,
         x2=>110,y2=>100,
         style=>{
            'stroke'=>"black"
         }
      );
         $path->path(
            d=>"M200 150 a150,30 0 1 0 0,70",
            style=>{stroke=>"purple",}
         );
         $path->path(
            d=>"M100 50 a1.2,5 0 1 1 20,0",
            style=>{stroke=>"purple",}
         );
         $path->path(
            d=>"M100 75 a1.2,10 0 1 1 30,0",
            style=>{stroke=>"purple",}
         );
         $path->path(
            d=>"M100 100 a1.2,2 0 1 1 10,0",
            style=>{stroke=>"purple",}
         );
      open OUT,">test.svg";
      print OUT $svg->xmlify;
      close OUT;
      system("/home/jfchen/159/FFproject/tools/draw/svg2xxx_release/svg2xxx -t pdf test.svg");

