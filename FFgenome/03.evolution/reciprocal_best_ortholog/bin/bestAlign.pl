#!/usr/bin/perl

=head1 Name

bestAlign.pl  --  choose the best hit for aligning

=head1 Description

This program can choose the best hit from aligning result,
accept psl m8 and solar format now. 
For psl, according to highest base matching rate
For m8, according to lowest e-value
For solar, according to highest aligning rate

=head1 Version

  Author: Fan Wei, fanw@genomics.org.cn
  Version: 1.0,  Date: 2007-2-7

=head1 Usage

  --fileformat   set input file format
  --cutoff       set a cut off to filter low quality alignments
  --verbose      output verbose information to screen  
  --help         output help information to screen  

=head1 Exmple

  perl bestAlign.pl example.psl -cutoff 0.5
  perl bestAlign.pl example.m8 -cutoff 1e-5
  perl bestAlign.pl example.solar -cutoff 0.5

=cut

use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname); 
use Data::Dumper;

my ($Fileformat,$Cutoff,$Verbose,$Help);
GetOptions(
	"fileformat:s"=>\$Fileformat,
	"cutoff:s"=>\$Cutoff,
	"verbose"=>\$Verbose,
	"help"=>\$Help
);
die `pod2text $0` if (@ARGV == 0 || $Help);

my %Data;

read_psl() if($Fileformat eq 'psl' || $ARGV[0] =~ /\.psl$/);
read_m8() if($Fileformat eq 'm8' || $ARGV[0] =~ /\.m8$/);
read_solar() if($Fileformat eq 'solar' || $ARGV[0] =~ /\.solar$/);

foreach my $qname (sort keys %Data) {
	print $Data{$qname}{line},"\n";
}

####################################################
################### Sub Routines ###################
####################################################


sub read_psl{
	while (<>) {
		chomp;
		my @temp = split /\t/;
		my $qname = $temp[9];
		my $value = $temp[0] / $temp[10]; ##精确匹配碱基率
		next if(defined $Cutoff && $value < $Cutoff);
		if (! exists $Data{$qname} || $Data{$qname}{value} < $value) {
			$Data{$qname}{value} = $value;
			$Data{$qname}{line} = $_;
		}
	}
}

sub read_m8{
	while (<>) {
		chomp;
		my @temp = split /\t/;
		my $qname = $temp[0];
		my $tname = $temp[1];
		next if($qname eq $tname);
		my $value = $temp[10]; ##E-value
		next if(defined $Cutoff && $value > $Cutoff);
		if (! exists $Data{$qname} || $Data{$qname}{value} > $value) {
			$Data{$qname}{value} = $value;
			$Data{$qname}{line} = $_;
		}
	}
}

sub read_solar{
	while (<>) {
		chomp;
		my @temp = split /\t/;
		my $qname = $temp[0];
		my $qlen = $temp[1];
		my $align_len;
		while ($temp[11] =~ /(\d+),(\d+);/g) {
			$align_len += $2 - $1 + 1;
		}
		my $value = $align_len / $qlen; ##align rate
		next if(defined $Cutoff && $value < $Cutoff);
		if (! exists $Data{$qname} || $Data{$qname}{value} < $value) {
			$Data{$qname}{value} = $value;
			$Data{$qname}{line} = $_;
		}
	}
}