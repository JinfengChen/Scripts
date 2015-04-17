#!/usr/bin/perl -w
#
########################################
#
# File Name:
#   sam2bed.pl
# 
# Description:
#   
# 
# Usage:
#   
# 
# Author:
#   Xi Wang, wang-xi05@mails.thu.edu.cn
# 
# Date:
#   Mon Nov  2 23:08:26 CST 2009
#
########################################

use strict;
my $usage = "$0 <.sam> <.bed> [.splice.bed]\n";
my $infile = shift || die $usage;
my $outfile = shift || die $usage;
my $splice = 0;
my $outfileSplice;
if ($outfileSplice = shift)
{
  $splice = 1;
}
open(IN, $infile) || die "Can't open $infile for reading!\n";
open(OUT, ">$outfile") || die "Can't open $outfile for writing!\n";
if ($splice) {
  open(OUTS, ">$outfileSplice") || die "Can't open $outfile for writing!\n";
}

my @col;
my ($chr, $s, $e, $str);
my $len;
my $mm;
my $i;
my $read;
while(<IN>)
{
    next if /^\@/;
    chomp;
    @col = split;
    $chr = $col[2];
    next if ($chr =~ /\*/);
    next if ($col[4] < 10); ### next if mapping quality < 10 to move low quality map
    $s = $col[3] - 1;
    $len = length($col[9]);
    $e = $s + $len;
    $str = $col[1] & 16;
    $str =~ s/0/+/;
    $str =~ s/16/-/;
    #/NM:i:(\d)/;
    /AS:i:(\d+)/;
    $mm = $1; ## alignment score
    $read=$col[0];
    if ($col[5] =~ /N/) # splice reads
    {
      if ($splice)
      {
        my ($lengths, $starts, $n);
        $starts = "0,";
        my $tmp;
        my @subcol = split /[NM]/, $col[5];
        $lengths = "$subcol[0],";
        $n = 1;
        $len = $subcol[0];
        for ($i=1; $i<@subcol; $i = $i + 2 )
        {
          $tmp = $len + $subcol[$i];
          $starts = "$starts$tmp,";
          if ($i + 1 > @subcol)
          {
            die "error read!\n";
          }
          $lengths = "$lengths$subcol[$i+1],";
          $n ++;
          $len = $len + $subcol[$i] + $subcol[$i+1];
          #print "$lengths\t$starts\n";
        }
        $e = $s + $len;
        print OUTS "$chr\t$s\t$e\tU$mm\t0\t$str\t-\t-\t-\t$n\t$lengths\t$starts\n";
      }
    }
    print OUT "$chr\t$s\t$e\t$read\t$mm\t$str\n";
}

close IN;
close OUT;
if ($splice)
{
  close OUTS;
}
