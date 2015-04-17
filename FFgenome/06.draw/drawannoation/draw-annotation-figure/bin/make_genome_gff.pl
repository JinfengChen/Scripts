#!/usr/bin/perl

=head1 Name

make_genome_gff.pl  --  make the genome gff file for glean

=head1 Description

The input is the fasta genome file, it invoke fastaDeal.pl 

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2006-12-6

=head1 Usage

  --verbose   output verbose information to screen  
  --help      output help information to screen  

=head1 Exmple

  perl ../bin/make_genome_gff.pl ./9311_Chr01.fa > ./9311_Chr01.fa.gff

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;

my ($Verbose,$Help);
GetOptions(
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV == 0 || $Help);

my $file = shift;

my $Result = `fastaDeal.pl $file -attr id:len`;

my @Res = split /\n/, $Result;

foreach  (@Res) {
	my ($id,$len) = ($1,$2) if(/(\S+)\s+(\S+)/);
	print "$id\tassembly\tscaffold\t1\t$len\t.\t+\t.\tScaffold $id\n"; 
}


####################################################
################### Sub Routines ###################
####################################################
