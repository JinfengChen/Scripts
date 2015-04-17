#!/usr/bin/perl
use strict;
use warnings;
use FindBin qw ($Bin);
use File::Basename qw(basename dirname); 
use Getopt::Long;

my %opt;
GetOptions(\%opt,"help:s");

if ($opt{help} or @ARGV < 1){
   print "In fasta sequence directory\n";
   print "Run :perl ../script_KaKs_calculator/runpipe.pl .";
}

my $qsub="/share/raid12/chenjinfeng/tools/bin/qsub-sge.pl";
my $pipe_kaks_file="/share/raid12/chenjinfeng/FFgenome/evolution/kaks_pairwise/bin/script_KaKs_calculator/pipe_kaks_file.pl";

my $dir=$ARGV[0];
my @files=glob("$dir/*.fa");

## write shell file 
my $kaksshell="pipe_kaks".".sh";
open OUT, ">$kaksshell" or die "can not open my shell out";
foreach (@files){
    print OUT "perl $pipe_kaks_file $_\n";
}
close OUT;

## run shell by qsub-sge.pl
`$qsub --resource vf=0.1G $kaksshell`;
my @inf=glob("$dir/*.inf");
foreach(@inf){
    `cat $_ >> $dir/KaKs.summary`;
    `rm $_`;
}


