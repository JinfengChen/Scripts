#!/usr/bin/perl

my @file=glob("*.blasttable");

foreach(@file){ 
    my $file=$_;
    if ($file=~/(.*)\.blasttable/){
       my $name=$1;
       $name=~tr/\_/\-/;
       system("/share/raid12/chenjinfeng/tools/inparanoid_4.1/table2inparanoid.pl $file > $name");
    }
}

