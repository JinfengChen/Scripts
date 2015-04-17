#!/usr/bin/perl

use Getopt::Long;

GetOptions(\%opt,"infile:s","outfile:s","help");

my $help=<<USAGE;
This program is design to check pfam domain in candidate sequences using hmmsearch (hmmer2.3).
For hmmer3.0, we can use --tblout or --domtblout to get a table format informations. 
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

my @pfam=glob("/share/raid12/chenjinfeng/FFgenome/evolution/Bigfamily/input/hmm2/*.hmm");


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
          #my $hmmrun="/share/raid1/genome/bin/hmmsearch --domE 1e-5 ".$pfam[$i]." "."$_.trans6";
          my $hmmrun="/share/raid1/genome/bin/hmmsearch ".$pfam[$i]." "."$_.trans6";
          my $str=`$hmmrun`;
          print "$str\n";
          #if ($str=~/\s---\n(.*)\n.*\n.*\n.*\n.*\n(.*)\n/){
          if ($str=~/\s---\n(.*)\n([\d\D]+)\s+-------\n([\d\D]+)\n\nAlignments of top-scoring domains/){
             #print "$3\n";
             my @temp_domains = split("\n", $3);
             foreach(@temp_domains){  
               #print "$_\n"; 
               my @temp_inf  = split(/\s+/, $_);
               if ($temp_inf[0]=~ /(.*)\_(\d+)/){
                   $flag=1;
                   $seqname=$1;
                   if ($2 > 3){
                      $strand="-";
                   }else{
                      $strand="+";
                   }
                   my $evalue=pop @temp_inf;
                   print OUT "$seqname\t$strand\t$domain\t$temp_inf[2]\t$temp_inf[3]\t$evalue\n";
               }
             }
          }
     }
}
close OUT;
system ("rm -R ./$filename.cut trans.log trans.log2");

