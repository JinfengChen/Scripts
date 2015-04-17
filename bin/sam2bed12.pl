#!/usr/bin/perl -w
#
########################################
#
# File Name:
#   sam2bed.pl
# 
# Description:
#   SAM format to BED format, splice junction reads considered.
# 
# Usage:
#   sam2bed.pl <.sam> <.bed>
# 
# Author:
#   Xi Wang, wang-xi05@mails.thu.edu.cn
# 
# Date:
#   Mon Nov  2 23:08:26 CST 2009
#
########################################

use strict;
my $usage = "$0 <.sam> <.bed>\n";
my $infile = shift || die $usage;
my $outfile = shift || die $usage;
open(IN, $infile) || die "Can't open $infile for reading!\n";
open(OUT, ">$outfile") || die "Can't open $outfile for writing!\n";

my (@col, @subcol);
my ($chr, $s, $e, $str);
my $len;
my $mm;
my $i;
my ($lengths, $starts, $n);
my $tmp;
my $junctionRead; # 1 if read map to junction database; 0 otherwise
my ($junc_s, $junc_l, $junc_r, $junc_e, $l_size, $r_size);

while(<IN>)
{
    next if /^\@/;
    chomp;
    @col = split /\t/;
    $chr = $col[2];
    next if ($chr =~ /\*/);
    if ($chr =~ /^(chr[0-9xyzmXYZM]+):(\d+)-(\d+)\|(\d+)-(\d+)$/)
    {
      $junctionRead = 1;
      $chr = $1;
      $junc_s = $2;
      $junc_l = $3;
      $junc_r = $4;
      $junc_e = $5;
      #print "$chr:$junc_s-$junc_l|$junc_r-$junc_e\n";
    }
    elsif ($chr =~ /^chr[0-9XYZMxyzm]+/)
    {
      $junctionRead = 0;
    }
    else
    {
      print STDERR "Can't recognize chr: $chr\n";
      next;
    }
    $s = $col[3] - 1;
    $str = $col[1] & 16;
    $str =~ s/0/+/;
    $str =~ s/16/-/;
    /NM:i:(\d)/;
    $mm = $1;

    $starts = "0,";
    @subcol = split /[NM]/, $col[5];
    if ($junctionRead == 1) # error: junction found in junction database
    {
      if(@subcol > 1)
      {
        print STDERR "Junction found in junction database.\n";
        next;
      }
      else
      {
        $n = 2;
        $s = $junc_s + $s;
        $e = $s + $junc_r - $junc_l + $subcol[0];
        $l_size = $junc_l - $s;
        $r_size = $e - $junc_r;
        $tmp = $junc_r - $s;
        $lengths = "$l_size,$r_size";
        $starts = "0,$tmp";
      }
    }
    else # junctionRead == 0
    {
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
    }
    print OUT "$chr\t$s\t$e\tU$mm\t0\t$str\t-\t-\t-\t$n\t$lengths\t$starts\n";
}

close IN;
close OUT;
