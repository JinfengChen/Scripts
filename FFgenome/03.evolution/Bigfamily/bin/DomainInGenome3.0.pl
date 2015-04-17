#!/usr/bin/perl

use Getopt::Long;

GetOptions(\%opt,"infile:s","outfile:s","help");

my $help=<<USAGE;
This program is design to check pfam domain in candidate sequences using hmmsearch (hmmer3.0).
perl DomainInGenome.pl -i ../input/super-scaffold -o super-scaffold.domain > log &
USAGE

if ($opt{help} or keys %opt < 1){
   print "$help\n";
   exit();
}

system ("/share/raid12/chenjinfeng/tools/bin/fastaDeal.pl -cuts 1 $opt{infile}");
my $filename=`basename $opt{infile}`;
chomp $filename;
my @file=glob("./$filename.cut/$filename.*");

my @pfam=glob("../input/hmm/*.hmm");


open OUT, ">$opt{outfile}" or die "$!";
foreach (@file){
     #print "$_\n";
     my %hash;
     my $flag;
     my $seqname;
     my $strand;
     my $tempout="$_".".trans6";
     system ("/share/raid1/genome/bin/transeq -sequence $_ -outseq $tempout -frame=6 > trans.log 2> trans.log2");
     for (my $i=0; $i<@pfam; $i++){
          my $domain=`basename $pfam[$i]`;
          chomp $domain;
          if ($domain=~/(.*)\.hmm/){
             $domain=$1;
          }
          
          my $hmmrun="/share/raid12/chenjinfeng/tools/hmmer3.0/bin/hmmsearch ".$pfam[$i]." "."$_.trans6";
          my $str=`$hmmrun`;
          print "$str\n";
          if ($str=~/\s-----------\n(.*)\n([\d\D]+)\s+----\n(.*)\n\n/){
          #if ($str=~/\s---\n(.*)\n([\d\D]+)\s+-------\n(.*)\n/){ hmmer2.3.2
               #print "$1\n";
               my @temp_plus = split(/\s+/, $1);
               my @temp_inf  = split(/\s+/, $3);
               #print "$temp_plus[9]\n";
               if ($temp_plus[9]=~ /(.*)\_(\d+)$/){
                   $flag=1;
                   $seqname=$1;
                   if ($2 > 3){
                      $strand="-";
                   }else{
                      $strand="+";
                   }
                   print OUT "$1\t$strand\t$domain\t$temp_inf[10]\t$temp_inf[11]\t$temp_plus[1]\n";
               }
          }
     }
}
close OUT;

system ("rm -R ./$filename.cut trans.log trans.log2");

